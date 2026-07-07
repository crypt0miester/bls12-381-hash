# bls381-hash-bench

Witness-assisted RFC 9380 hash-to-curve for BLS12-381 in a Solana SBF program,
with a mollusk compute-unit benchmark. Both the G1 (min-sig) and G2 (min-pk)
hashes run using only syscalls active on Solana mainnet today. no `big_mod_exp`
(SIMD-0529) and no map-to-curve syscall.

## Approach

Everything `big_mod_exp` would compute (field inverses, Legendre symbols, square
roots) is expensive to compute but cheap to verify, so the caller supplies each
as instruction data (a witness) and the program checks it with one or two field
multiplications:

- sqrt: witness `y`, check `y^2 == gx`
- inverse: witness `t`, check `v*t == 1`

A wrong witness aborts the instruction, and no witness can steer the output, so
the hash stays a pure function of the message. Two reductions on top of the basic
scheme:

- one sqrt witness per SSWU branch:  `g(x2) = (Z u^2)^3 g(x1)` with `Z` a
  non-residue, so `gx1`/`gx2` always have opposite quadratic character and a
  single sqrt proves its own branch (the `gx1 == 0` degenerate is rejected on the
  x2 branch to match blst)
- Montgomery pair inversion:  one witness `w = (a*b)^-1` pins both inverses of a
  same-stage pair via a single product check

## Layout

- `program/src/g1_msig.rs`:  hash-to-G1, min-sig suite (`BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_`)
- `program/src/g2_msig.rs`:  hash-to-G2, min-pk suite (G2 equivalent)
- `program/src/{g1,g2}_consts.rs`:  SSWU / isogeny / psi constants, extracted from zkcrypto `bls12_381`
- `program/src/lib.rs`:  instruction dispatch (tagged probes + the hash pipelines)
- `bench/tests/bench.rs`:  mollusk benchmark, blst cross-checks, corrupt-witness rejection

## Build and run

```
cd program && cargo build-sbf --arch v3
cd ../bench && cargo test -- --nocapture
```

Requires the Solana platform tools (`cargo build-sbf`). The tests assert
byte-equality with `blst` on the standard POP suites for every stage and check
that a corrupted witness aborts rather than producing a different point.

## Measured (mollusk 0.13.4, agave 4.0, SBF v3)

```
hash_to_G1 (min-sig):  ~363k CU   338 B witness
hash_to_G2 (min-pk):   ~588k CU   482 B witness
```

## Status

Experimental. The witnessed hash is novel enough to warrant 
a hostile review before it is used anywhere consensus depends on it.
