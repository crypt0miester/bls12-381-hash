# bls12-381-hash

Witness-assisted RFC 9380 hash-to-curve for BLS12-381, for Solana SBF programs.
`no_std`, allocation-free on-chain, using only syscalls active on mainnet today,
no `big_mod_exp` (SIMD-0529) and no map-to-curve syscall. Everything expensive
(inverses, Legendre symbols, square roots) rides in as witness data and is
verified with a multiplication or two; the byte-minimal G2 layout shrinks the
witness to 97 bytes by shipping half of each square root and recomputing the
one remaining inverse in-program. Host-side witness generation ships in the
same crate.

```rust
use bls381_hash::{dst, hash_to_g1, hash_to_g2_compact_xgcd};

// on-chain: DST is a runtime parameter, the payload is the witness bytes
// followed by the message
let point = hash_to_g1(dst::G1_RO, payload)?; // Vec<u8>, the 96-byte G1 point

// the byte-minimal G2 layout: a 97-byte blob (branch flags plus the two
// square-root halves), the batched inverse recomputed in-program
let point = hash_to_g2_compact_xgcd(dst::G2_RO, payload)?;
```

```rust
// off-chain (host): build the witness for a message
let witness = bls381_hash::witness::g1::generate(message);
let witness = bls381_hash::witness::g2::generate_compact_xgcd(message);
```

## Features

| feature | pulls in |
|---|---|
| `ro` (default) | `g1-ro` + `g2-ro`, the blst-compatible pair |
| `g1-ro`, `g2-ro` | standard `_SSWU_RO_POP_` pipelines; `g2-ro` includes the compact and xgcd blob layouts |
| `g1-nu`, `g2-nu` | RFC 9380 encode_to_curve variants; not random-oracle suites, see the NU note in BENCHMARKS.md |
| `modexp` | big_mod_exp-assisted G1 path, needs SIMD-0529 |
| `wide-witness` | 674 B G2 blob, ~14k CU cheaper; for 4 KiB (SIMD-0296) transactions |
| `full` | everything above |

The `lib/` crate (`bls381-hash`) is the product; `program/` is an SBF
tag-dispatch fixture and `bench/` the mollusk benchmark.

## Cost

Measured with mollusk 0.13.4 on the agave 4.0 stack, SBF v3:

| pipeline | CU | witness | compatibility |
|---|---|---|---|
| hash_to_G1 (RO, min-sig) | ~129k | 338 B | `_SSWU_RO_POP_`, byte-equal to blst |
| hash_to_G2 (RO, min-pk) | ~253k | 578 B | `_SSWU_RO_POP_`, byte-equal to blst |
| hash_to_G2 (RO, compact) | ~470k | 145 B | same suite, same output bytes |
| hash_to_G2 (RO, xgcd) | ~583k | 97 B | same suite, same output bytes |
| hash_to_G2 (RO, `wide-witness`) | ~241k | 674 B | same suite, bigger blob |

An end-to-end min-pk BLS verify (hash_to_G2 plus the pairing syscall) lands
around 308k CU with the default blob, 525k compact, 638k xgcd; for scale, a
naive port of zkcrypto `bls12_381` costs 11.3M CU for G1 and 46.5M CU for
G2. The NU variants, stage costs, syscall pricing, and the optimization
notes live in [BENCHMARKS.md](BENCHMARKS.md).

## Picking a blob

The G2 layouts trade witness bytes for CU, and the trade only ever binds on
bytes: a legacy transaction caps at 1,232 of them while the CU budget holds
1.4M. The 97 B xgcd blob is the chosen default: it is small enough that
every BLS-verifying instruction surveyed (committee votes, certifications
with a 648 B storage proof, key registration with a second signer) carries
its witness inline in one legacy transaction, which retires stage-then-
consume flows entirely. The 145 B blob buys ~112k CU back for 48 bytes when
the transaction has room; the 578 B blob is the CU floor for consumers with
no byte pressure. All three run the same checks and produce byte-identical
output, so the choice is per call site. Once SIMD-0529 (`big_mod_exp`)
activates, the square-root halves become computable on-chain too and the
witness drops to zero bytes.

## Approach

Field inverses, Legendre symbols, and square roots are expensive to compute
but cheap to verify, so the caller supplies each as instruction data and the
program checks it:

- sqrt: witness `y`, check `y^2 == gx`
- inverse: witness `t` in Montgomery form, check `v*t == 1`
- quotient: witness the result `q` of `n/d` directly, check `q*d == n`

The compact layouts restructure the same checks to cut bytes:

- a square root travels as its real part alone: with `gx = g0 + g1 i`, the
  imaginary half of `y^2 == gx` forces `c1 = g1/(2 c0)`, and the real half,
  checked as `(c0 + c1)(c0 - c1) == g0`, pins the root
- every inverse the pipeline needs hides behind one 48-byte batched witness
  `w = (e1...e8)^-1`, pinned by a single product check; Fp2 inverses ride
  their norms, and a zero anywhere fails the product
- the xgcd layout drops even that witness: inversion is gcd-shaped, so the
  program recomputes `w` with a binary extended gcd, shift-and-subtract
  only, no multiplies (~112k CU)
- the square-root halves stay bytes because sqrt mod p is
  exponentiation-shaped, with no multiply-free algorithm known

A wrong witness aborts the instruction, and no witness can steer the output,
so the hash stays a pure function of the message. How each cost was carved
down, what was measured and abandoned, and which knobs remain are in
[BENCHMARKS.md](BENCHMARKS.md).

## Layout

- `lib/src/fp.rs`, `fp2.rs`: Fp and Fp2 arithmetic (ps30 product-scanning Montgomery; the straight-line bodies come from `tools/gen_ps30.py`, which also proves the column bounds and simulates the emitted algorithm against a reference reduction), plus the binary-xgcd inverse behind the witness-free batch
- `lib/src/g1.rs`, `g2.rs`: RO and NU pipelines, the compact and xgcd G2 layouts, cofactor clearing, host witness generation
- `lib/src/consts_g1.rs`, `consts_g2.rs`: SSWU, isogeny, psi, and adapted constants, with the compile-time 2^390-domain lift
- `lib/src/lib.rs`: public API, feature gates, `dst` module
- `program/src/lib.rs`: SBF tag-dispatch fixture
- `bench/tests/`: mollusk benchmarks, blst cross-checks, soundness tests
- `BENCHMARKS.md`: stage costs, syscall pricing, optimization log and open knobs

## Build and run

```
cd program && cargo build-sbf --arch v3
cd ../bench && cargo test -- --nocapture
```

Requires the Solana platform tools. The standard suites assert byte-equality
with blst at every stage. A corrupted witness must abort, and supplying the
other square root must not change the output point.

## Status

Experimental. The witnessed hash is novel enough to warrant a hostile review
before it is used anywhere consensus depends on it.

## License

MIT. The SSWU, isogeny, and psi constants in `lib/src/consts_g{1,2}.rs` were
extracted from zkcrypto [`bls12_381`](https://github.com/zkcrypto/bls12_381)
(MIT/Apache-2.0); the map constructions follow Wahby-Boneh, eprint 2019/403.
