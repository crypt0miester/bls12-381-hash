# bls12-381-hash

Witness-assisted RFC 9380 hash-to-curve for BLS12-381, for Solana SBF programs.
`no_std`, allocation-free on-chain, using only syscalls active on mainnet today,
no `big_mod_exp` (SIMD-0529) and no map-to-curve syscall. Everything expensive
(inverses, Legendre symbols, square roots) rides in as witness data and is
verified with a multiplication or two; the byte-minimal G2 layout shrinks the
witness to 96 bytes by shipping half of each square root, riding the branch
bits on their parity, and recomputing the one remaining inverse in-program.
Host-side witness generation ships in the same crate.

```rust
use bls381_hash::{dst, hash_to_g1, hash_to_g2_compact_parity};

// on-chain: DST is a runtime parameter, the payload is the witness bytes
// followed by the message
let point = hash_to_g1(dst::G1_RO, payload)?; // Vec<u8>, the 96-byte G1 point

// the byte-minimal G2 layout: a 96-byte blob (the two square-root halves,
// each branch flag riding its half's parity), the batched inverse
// recomputed in-program by divsteps
let point = hash_to_g2_compact_parity(dst::G2_RO, payload)?;
```

```rust
// off-chain (host): build the witness for a message
let witness = bls381_hash::witness::g1::generate(message);
let witness = bls381_hash::witness::g2::generate_compact_parity(message);
```

## Features

| feature | pulls in |
|---|---|
| `ro` (default) | `g1-ro` + `g2-ro`, the blst-compatible pair |
| `g1-ro`, `g2-ro` | standard `_SSWU_RO_POP_` pipelines; `g2-ro` includes the compact, xgcd and parity blob layouts |
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
| hash_to_G2 (RO, xgcd) | ~509k | 97 B | same suite, same output bytes |
| hash_to_G2 (RO, parity) | ~509k | 96 B | same suite, same output bytes |
| hash_to_G2 (RO, `wide-witness`) | ~241k | 674 B | same suite, bigger blob |

An end-to-end min-pk BLS verify (hash_to_G2 plus the pairing syscall) lands
around 308k CU with the default blob, 525k compact, 564k parity; for scale,
a naive port of zkcrypto `bls12_381` costs 11.3M CU for G1 and 46.5M CU for
G2. The NU variants, stage costs, syscall pricing, and the optimization
notes live in [BENCHMARKS.md](BENCHMARKS.md).

## Picking a blob

The G2 layouts trade witness bytes for CU, and the trade only ever binds on
bytes: a legacy transaction caps at 1,232 of them while the CU budget holds
1.4M. The 96 B parity blob is the chosen default: it is small enough that
every BLS-verifying instruction surveyed (committee votes, certifications
with a 648 B storage proof, key registration with a second signer) carries
its witness inline in one legacy transaction, which retires stage-then-
consume flows entirely, and it is the floor of this witness family (the
two root halves are pure computational advice for the square roots, and
sqrt mod p has no multiply-free algorithm). The 145 B blob buys ~39k CU
back for 49 bytes when the transaction has room; the 578 B blob is the CU
floor for consumers with no byte pressure. All of them run the same checks
and produce byte-identical output, so the choice is per call site. Once
SIMD-0529 (`big_mod_exp`) activates, the square-root halves become
computable on-chain too and the witness drops to zero bytes.

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
  program recomputes `w` with Bernstein-Yang divsteps over 30-bit lanes
  (~39k CU; `tools/check_divsteps.py` proves the batching against the
  paper's divstep and checks the code lane for lane against Fermat)
- the parity layout drops the flags byte too: the verifier re-canonicalizes
  the root sign anyway, so each branch bit rides its root half's parity,
  and the witness becomes unique (the other root reads as a branch lie and
  aborts)
- the square-root halves stay bytes because sqrt mod p is
  exponentiation-shaped, with no multiply-free algorithm known

A wrong witness aborts the instruction, and no witness can steer the output,
so the hash stays a pure function of the message. How each cost was carved
down, what was measured and abandoned, and which knobs remain are in
[BENCHMARKS.md](BENCHMARKS.md).

## Algorithms and references (26 total)

Curve and hash construction (15):

- RFC 9380 hash-to-curve: `expand_message_xmd`, `hash_to_field`, `sgn0`,
  the RO/NU suites and DSTs, the isogeny map constants (appendix E), and
  the square-root rules for p = 3 mod 4 (appendix I).
  [rfc-editor.org/rfc/rfc9380](https://www.rfc-editor.org/rfc/rfc9380)
- The map to the curve is simplified SWU: Shallue and van de Woestijne,
  ANTS 2006 ([doi 10.1007/11792086_36](https://doi.org/10.1007/11792086_36));
  Ulas ([arXiv:0706.1448](https://arxiv.org/abs/0706.1448)); simplified by
  Brier, Coron, Icart, Madore, Randriam and Tibouchi
  ([eprint 2009/340](https://eprint.iacr.org/2009/340)).
- Wahby and Boneh: SSWU on an isogenous curve for A B = 0, the 11- and
  3-isogeny construction, the opposite quadratic character of gx1 and gx2
  that the one-root branch proof leans on, and the effective-cofactor G1
  clearing. [eprint 2019/403](https://eprint.iacr.org/2019/403)
- G2 cofactor clearing: Budroni and Pintore
  ([eprint 2017/419](https://eprint.iacr.org/2017/419)), the endpoint of
  the psi-based line through Scott, Benger, Charlemagne, Dominguez Perez
  and Kachisa ([eprint 2008/530](https://eprint.iacr.org/2008/530)) and
  Fuentes-Castaneda, Knapp and Rodriguez-Henriquez, SAC 2011
  ([doi 10.1007/978-3-642-28496-0_25](https://doi.org/10.1007/978-3-642-28496-0_25)),
  built on the psi (untwist-Frobenius) endomorphism of Galbraith and Scott
  ([eprint 2008/117](https://eprint.iacr.org/2008/117)). RFC 9380 froze
  this map for the suite, so the clearing is normative, not an
  optimization slot.
- The curve family: Barreto, Lynn and Scott
  ([eprint 2002/088](https://eprint.iacr.org/2002/088)); BLS12-381 itself:
  [Bowe](https://electriccoin.co/blog/new-snark-curve/); parameters and the
  Zcash point serialization:
  [draft-irtf-cfrg-pairing-friendly-curves](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/).
- The min-pk verify fixture: Boneh, Lynn and Shacham signatures
  ([doi 10.1007/s00145-004-0314-9](https://doi.org/10.1007/s00145-004-0314-9))
  per
  [draft-irtf-cfrg-bls-signature](https://datatracker.ietf.org/doc/draft-irtf-cfrg-bls-signature/).

Field arithmetic (11):

- Montgomery multiplication and reduction: Montgomery, "Modular
  multiplication without trial division", Math. Comp. 44 (1985)
  ([doi 10.1090/S0025-5718-1985-0777282-X](https://doi.org/10.1090/S0025-5718-1985-0777282-X)).
- The word-level Montgomery variants (CIOS and friends, the measured
  baseline ps30 beats): Koc, Acar and Kaliski, "Analyzing and comparing
  Montgomery multiplication algorithms", IEEE Micro 16 (1996)
  ([pdf](https://cetinkayakoc.net/docs/j37.pdf)).
- Product-scanning (column) ordering: Comba, IBM Systems Journal 29 (1990)
  ([doi 10.1147/sj.294.0526](https://doi.org/10.1147/sj.294.0526)).
- Unsaturated limbs (lanes narrower than the word, carries deferred to
  provable column bounds): the school of Bernstein's Curve25519 radix-25.5
  arithmetic, PKC 2006
  ([pdf](https://cr.yp.to/ecdh/curve25519-20060209.pdf)). The 30-bit-lane
  ps30 instantiation and its column-bound proofs are this crate's
  (`tools/gen_ps30.py`).
- Fp2 multiplication is Karatsuba
  ([Karatsuba and Ofman, 1962](https://en.wikipedia.org/wiki/Karatsuba_algorithm)),
  squaring the (a + b)(a - b) complex form.
- The isogeny polynomials run with preprocessed (adapted) coefficients:
  Knuth, The Art of Computer Programming, vol. 2, section 4.6.4.
- The one shared inverse behind every witness batch: Montgomery
  simultaneous inversion, "Speeding the Pollard and elliptic curve methods
  of factorization", Math. Comp. 48 (1987)
  ([doi 10.1090/S0025-5718-1987-0866113-7](https://doi.org/10.1090/S0025-5718-1987-0866113-7)).
- The in-program inverse: Bernstein and Yang divsteps ("safegcd") and its
  iteration bound ([eprint 2019/266](https://eprint.iacr.org/2019/266));
  the 30-bit lane batching and the f(f^2 - 2) = -1/f mod 2^6 identity
  follow the
  [libsecp256k1 implementation notes](https://github.com/bitcoin-core/secp256k1/blob/master/doc/safegcd_implementation.md);
  the port is mirrored and bound-checked by `tools/check_divsteps.py`.
- Host witness generation: Fermat inverse and the Euler criterion, and
  fixed 4-bit-window exponentiation, per the
  [Handbook of Applied Cryptography](https://cacr.uwaterloo.ca/hac/)
  (chapters 2 and 14); Fp2 square roots by the complex (norm) method: Adj
  and Rodriguez-Henriquez
  ([eprint 2012/685](https://eprint.iacr.org/2012/685)).

The witness-verified pipeline itself, the half-root transport
(`c1 = g1/(2 c0)` with the difference-of-squares check), the single
homogenized inverse batch that closes before any inverse exists, and the
parity-encoded branch bits have no external reference; they are this
crate's constructions, argued in [BENCHMARKS.md](BENCHMARKS.md).

## Layout

- `lib/src/fp.rs`, `fp2.rs`: Fp and Fp2 arithmetic (ps30 product-scanning Montgomery; the straight-line bodies come from `tools/gen_ps30.py`, which also proves the column bounds and simulates the emitted algorithm against a reference reduction), plus the divsteps inverse behind the witness-free batch (mirrored and bound-checked by `tools/check_divsteps.py`)
- `lib/src/g1.rs`, `g2.rs`: RO and NU pipelines, the compact, xgcd and parity G2 layouts, cofactor clearing, host witness generation
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
other square root must not change the output point (in the parity layout it
must abort outright: the root's parity is the branch bit).

## Status

Experimental. The witnessed hash is novel enough to warrant a hostile review
before it is used anywhere consensus depends on it.

## License

MIT. The SSWU, isogeny, and psi constants in `lib/src/consts_g{1,2}.rs` were
extracted from zkcrypto [`bls12_381`](https://github.com/zkcrypto/bls12_381)
(MIT/Apache-2.0); the map constructions follow Wahby-Boneh, eprint 2019/403.
