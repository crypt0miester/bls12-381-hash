# bls381-hash-bench

RFC 9380 hash-to-curve for BLS12-381 inside a Solana SBF program, using only
syscalls that are active on mainnet today. No `big_mod_exp` (SIMD-0529) and no
map-to-curve syscall. Everything expensive rides in as witness data and gets
verified with a multiplication or two. Ships with a mollusk compute-unit
benchmark cross-checked against blst.

## Results

Measured with mollusk 0.13.4 on the agave 4.0 stack, SBF v3.

| pipeline | CU | witness | compatibility |
|---|---|---|---|
| hash_to_G1 (RO, min-sig) | ~245k | 338 B | `_SSWU_RO_POP_`, byte-equal to blst |
| hash_to_G2 (RO, min-pk) | ~484k | 482 B | `_SSWU_RO_POP_`, byte-equal to blst |
| encode_to_G1 (NU) | ~183k | 193 B | `_SSWU_NU_POP_`, byte-equal to blst encode |
| encode_to_G2 (NU) | ~326k | 289 B | `_SSWU_NU_POP_`, byte-equal to blst encode |
| hash_to_G1 (SvdW) | ~162k | 242 to 434 B | custom suite |
| hash_to_G2 (SvdW) | ~426k | 482 to 866 B | custom suite |

For scale, a naive port of zkcrypto `bls12_381` costs 11.3M CU for G1 and
46.5M CU for G2, and a single 381-bit field multiplication bottoms out around
3.4k CU on sbpf.

The NU suites hash with a single map (RFC 9380 encode_to_curve). Note that
the CFRG BLS signature draft registers only hash_to_curve (RO) ciphersuites,
and RFC 9380 limits encode_to_curve to applications whose security analysis
does not rely on a random oracle. Using NU for BLS rests on the argument in
section 5 of Wahby-Boneh ([eprint 2019/403](https://eprint.iacr.org/2019/403))
and the BCI+10 reference there: hashing onto a constant fraction of the group
suffices for unforgeability. That makes it a deliberate protocol choice, and
the hash must never be reused for anything that actually needs a random
oracle.

## Approach

Field inverses, Legendre symbols, and square roots are expensive to compute
but cheap to verify, so the caller supplies each as instruction data and the
program checks it:

- sqrt: witness `y`, check `y^2 == gx`
- inverse: witness `t` in Montgomery form, check `v*t == 1`
- non-square (SvdW branches): witness `s` with `s^2 = xi * f(x)` for a fixed
  non-residue `xi`

A wrong witness aborts the instruction, and no witness can steer the output,
so the hash stays a pure function of the message.

## Optimizations that landed

- One sqrt witness per SSWU branch. `g(x2) = (Z u^2)^3 g(x1)` with `Z` a
  non-residue (section 4.1 of the paper), so a single root proves its own
  branch.
- Montgomery pair inversion. One witness `w = (a*b)^-1` pins both inverses of
  a same-stage pair.
- Knuth-adapted polynomial evaluation (TAOCP 4.6.4 preprocessing). The iso-11
  runs in 27 multiplications instead of 51, the iso-3 in 5 plus a square
  instead of 11. Constants were derived and expansion-checked offline.
- Small constants (`Z = 11`, `Z2 = -(2+i)`, `A' = 240i`, the `2^256` fold)
  are multiplied with addition chains instead of field muls.

## SvdW variant (custom suite)

A direct Shallue-van de Woestijne map onto the curve (section 3 of the paper,
`u0 = -3` on E1 and `u0 = -1` on E2), which skips the isogeny entirely.
Branch selection is proven rather than trusted: SvdW takes the smallest j
with `f(x_j)` square, so claiming branch j requires a non-squareness proof
for each earlier branch plus the sqrt witness for its own. The output matches
no standard suite, and the NU suites beat it on G2 anyway, so this one stays
a finding.

## Measured dead ends

All of these are kept runnable (`mul_bench` variants and
`bench/tests/modexp_cost.rs`):

- A dedicated 32-bit squaring costs 4.4k CU against 3.4k for the general
  multiply. The doubling pass and carry ripples outweigh the saved cross
  products.
- A lazy-carry 30-bit comba multiply costs 7.1k CU. Bounds checks and limb
  conversion bury the deferred carries.
- A fully unrolled cios32 lands at 3.38k CU, identical to the loop, because
  LLVM already unrolls it. Around 3.4k really is the ISA floor; we tried
  seven implementations.
- Multiplication through modexp identities fails because `big_mod_exp` prices
  flat at about 1.7k CU regardless of exponent, and the Montgomery correction
  eats the win.
- The paper's fused sqrt-ratio and projective evaluation (section 4.2) trade
  cheap witnessed exponentiations for 3.4k muls. The CPU cost model that the
  hash-to-curve literature optimizes for is inverted on sbpf.
- Delayed Montgomery reduction across Fp2 (`mul2_lazy`): three wide products
  and two reductions instead of three of each, 17% fewer multiply-accumulates,
  and still 0.5% slower on the full G2 pipeline. The separate wide passes cost
  more than the saved reduction. This also settles base-field Karatsuba, which
  restructures the same way for a smaller saving.

## Further optimizations

Open knobs, in rough order of interest:

- The modexp path (tags 30 to 33) runs hash_to_G1 in ~295k CU with zero
  witness bytes, against ~245k plus 338 B for the witnessed path. It needs
  `big_mod_exp` (SIMD-0529, merged but not active). Once 0529 activates, a
  transaction that is byte-bound rather than CU-bound should prefer it.
- The min-pk verify transaction is closer to byte-bound than CU-bound (537k
  of the 1.4M CU ceiling, but witness plus keys eat real transaction space).
  Batching all five G2 inverses behind one witness would save 192 B for
  roughly 17k CU. Not implemented; the right trade depends on the consumer.
- G2 cofactor clearing costs ~45k CU across ~140 g2 add syscalls. The
  Budroni-Pintore chain is the best known construction; the cost is syscall
  pricing, not structure. The verify path feeds the hash into the pairing
  uncompressed, so no decompression cost hides there.

## Layout

- `program/src/g1_msig.rs`: G1 RO and NU pipelines, Fp arithmetic
- `program/src/g2_msig.rs`: G2 RO and NU pipelines, Fp2 arithmetic, cofactor clearing
- `program/src/{g1,g2}_svdw.rs`: SvdW variants
- `program/src/{g1,g2}_consts.rs`: SSWU, isogeny, psi, and adapted constants
- `program/src/lib.rs`: instruction dispatch
- `bench/tests/`: mollusk benchmarks, blst cross-checks, corrupt-witness rejection

## Build and run

```
cd program && cargo build-sbf --arch v3
cd ../bench && cargo test -- --nocapture
```

Requires the Solana platform tools. The standard suites assert byte-equality
with blst at every stage. SvdW checks against a host-side reference. A
corrupted witness must abort, and supplying the other square root must not
change the output point.

## Status

Experimental. The witnessed hash is novel enough to warrant a hostile review
before it is used anywhere consensus depends on it.

## License

MIT. The SSWU, isogeny, and psi constants in `{g1,g2}_consts.rs` were
extracted from zkcrypto [`bls12_381`](https://github.com/zkcrypto/bls12_381)
(MIT/Apache-2.0); the map constructions follow Wahby-Boneh, eprint 2019/403.
