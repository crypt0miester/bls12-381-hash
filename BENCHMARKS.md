# Benchmark details

Measured with mollusk 0.13.4 on the agave 4.0 stack, SBF v3, against the
`program/` fixture; the bench suite asserts byte-equality with blst at every
stage.

Known issue (2026-07-13): the naive zkcrypto reference rows (tags 0/1)
produce wrong points on the current platform-tools; both naive tests fail,
silently and at plausible CU, and their historical numbers (11.3M G1 /
46.5M G2) date from the earlier toolchain. Suspected u128 miscompilation;
the crate's own paths use no u128 anywhere, are unaffected, and stay blst
byte-verified.

## Pipelines

| pipeline | CU | witness | compatibility |
|---|---|---|---|
| hash_to_G1 (RO, min-sig) | ~129k | 338 B | `_SSWU_RO_POP_`, byte-equal to blst |
| hash_to_G2 (RO, min-pk) | ~253k | 578 B | `_SSWU_RO_POP_`, byte-equal to blst |
| hash_to_G2 (RO, compact) | ~470k | 145 B | same suite, same output bytes |
| hash_to_G2 (RO, xgcd) | ~583k | 97 B | same suite, same output bytes |
| hash_to_G2 (RO, `wide-witness`) | ~241k | 674 B | same suite, bigger blob |
| encode_to_G1 (NU) | ~99k | 193 B | `_SSWU_NU_POP_`, byte-equal to blst encode |
| encode_to_G2 (NU) | ~173k | 385 B | `_SSWU_NU_POP_`, byte-equal to blst encode |

The default and `wide-witness` G2 rows are the same pipeline with a
different tv2 witness layout: the default 578 B blob pins both tv2 inverses
behind one pair inversion witness, while `wide-witness` ships both inverses
directly (96 more bytes, ~12k CU less) for when SIMD-0296 4 KiB
transactions land. The compact and xgcd rows are the byte-bound layouts,
detailed below. An end-to-end min-pk BLS verify (hash_to_G2 plus the
pairing syscall) lands around 308k CU (297k wide, 525k compact, 638k xgcd).
A 381-bit field multiplication costs ~1.9k CU with the ps30
product-scanning multiplier (the textbook 32-bit CIOS form bottoms out
around 3.3k).

The NU suites hash with a single map (RFC 9380 encode_to_curve). Note that
the CFRG BLS signature draft registers only hash_to_curve (RO) ciphersuites,
and RFC 9380 limits encode_to_curve to applications whose security analysis
does not rely on a random oracle. Using NU for BLS rests on the argument in
section 5 of Wahby-Boneh ([eprint 2019/403](https://eprint.iacr.org/2019/403))
and the BCI+10 reference there: hashing onto a constant fraction of the group
suffices for unforgeability. That makes it a deliberate protocol choice, and
the hash must never be reused for anything that actually needs a random
oracle.

## Compact witness layouts

Two G2 layouts restructure the same checks to shrink the blob for
byte-bound transactions (hash tags 48/49, end-to-end verify tags 54/55,
exact numbers at bench commit of 2026-07-13):

| layout | blob | hash CU | e2e verify k=14 / k=20 |
|---|---|---|---|
| default | 578 B | 252,916 | 308,209 / 294,568 |
| compact | 145 B | 470,108 | 525,373 / 511,756 |
| xgcd | 97 B | 582,583 | 637,851 / 624,234 |

The compact blob is one flags byte (the two SSWU branch bits), the real
halves of the two square roots, and one batched inverse witness. A root
travels as c0 alone: the imaginary half of `y^2 == gx` forces
`c1 = g1/(2 c0)`, the real half, checked as `(c0 + c1)(c0 - c1) == g0`,
pins it, and both roots collapse to one output through the sgn0 rule (the
sign read is c0's parity, since a zero c0 cannot pass the batch). Every
inverse the pipeline needs, the 2c0 divisors, the tv2 pair, and the
per-map iso-3 denominators, hides behind `w = (e1..e8)^-1` with one
product check; Fp2 inverses ride their norms (`norm z = 0` only at
`z = 0` since -1 is a non-residue), so a zero anywhere fails the product.
The iso-3 denominators join the batch before any inverse exists by
evaluating them homogenized over the fractional candidate `x' = n/tv2`
(`compact_map_parts`, shared verbatim by the verifier and the host
generator so the batch layout is defined once), each map lands on E on
its own, and the E-side addition rides the g2 add syscall instead of a
slope witness.

The xgcd layout drops the last witness: inversion is gcd-shaped, so the
batch product is inverted in-program by binary extended gcd
(`fp.rs::inv_xgcd`, shift-and-subtract only, no multiplies, a 768-pass
cap; ~112k CU, mildly message-dependent). Square roots have no
multiply-free algorithm, so the two c0 halves are the transported floor,
97 B with the flags byte. A branch lie either desyncs the batch from `w`
or leaves the real-part check unsatisfiable (the wrong branch's gx is a
non-square); the soundness sweeps cover both layouts bit by bit,
including batch-consistent lies (`generate_compact_*_steered`) and the
other square root (`flip_first_root`, which also negates the trailing
inverse when the layout carries one).

Two measure-zero deviations beyond the pipeline's existing aborts: a gx
that is zero on the x1 branch, and a root with zero real part, both abort
instead of hashing (probability ~2^-380, unreachable for SHA-derived
inputs).

## Stage costs

Tags 46/47, cumulative prefixes of the witnessed pipelines:

| stage | min-sig G1 | min-pk G2 |
|---|---|---|
| hash_to_field | 9.7k | 19.7k |
| both SSWU maps | 34.4k | 100.1k |
| E' add + isogeny | 69.7k | 78.0k |
| clear_cofactor + validate | 14.9k | 55.1k |

The isogeny evaluation dominates the field work (soundness requires
evaluating all four polynomials at the summed point, so no witness can
shortcut it), and the clearing stages are syscall pricing.

## Optimizations that landed

- One sqrt witness per SSWU branch. `g(x2) = (Z u^2)^3 g(x1)` with `Z` a
  non-residue (section 4.1 of the Wahby-Boneh paper), so a single root
  proves its own branch.
- Montgomery pair inversion. One witness `w = (a*b)^-1` pins both inverses of
  a same-stage pair.
- Knuth-adapted polynomial evaluation (TAOCP 4.6.4 preprocessing). The iso-11
  runs in 27 multiplications instead of 51, the iso-3 in 5 plus a square
  instead of 11. Constants were derived and expansion-checked offline.
- Small constants (`Z = 11`, `Z2 = -(2+i)`, `A' = 240i`, the `2^256` fold)
  are multiplied with addition chains instead of field muls.
- Bare Montgomery reduction out of Montgomery form. `from_mont` skips the
  product loop of a multiply by one, about half the multiply-accumulates.
- Unreduced sums: the mont_mul operand bound is 2^384, not p, so the Fp2
  Karatsuba sums and the iso polynomial stage sums skip their conditional
  subtraction (an add_mod measures ~120 CU, its unreduced form about 60;
  bounds at the stage macro and iso evaluations).
- Syscall result buffers are `MaybeUninit`: the pointer escapes into the
  syscall, so LLVM cannot drop the zero-init on its own, and cofactor
  clearing makes ~70 (G1) / ~145 (G2) such calls per hash (about 1.4k and
  2.8k CU respectively).
- Product-scanning Montgomery multiplication over 30-bit lanes ("ps30",
  `fp.rs`): columns accumulate into a single u64 with no per-product masking,
  carrying or stores, and the modulus lanes ride as immediate operands. One
  multiply costs ~1.93k CU against ~3.3k for 32-bit CIOS. The radix moves to
  `R = 2^390 mod p`; every Montgomery constant is lifted at compile time
  (`consts_g{1,2}.rs`). SBPF v3 has no wide multiply (the v2-only PQR class
  with UHMUL was dropped), so a u128 product lowers to a ~70 CU `__multi3`
  call: 64-bit limbs are off the table. Thirteen lanes is the floor for this
  ISA: 12x32 breaks the a*b < pR bound (and 32-bit masks are not valid
  immediates), and 31-bit lanes keep the same product count as 30 while
  overflowing far more columns. The earlier 29-bit form (ps29, 14 lanes,
  2.17k CU) kept every column sum below 2^64 for free; ps30 spends one lane
  less (392 -> 338 products) and pays with a banked-carry spill in the three
  tallest columns, whose exact bounds `tools/gen_ps30.py` tracks and proves.
- Dedicated ps30 squaring (doubled cross products) for the Fp squarings in
  the SSWU map and iso-11 chain.
- Witnesses answer exactly one question each. The E' addition slope is
  witnessed directly (`lambda * dx == dy`) and the iso-map output
  coordinates are pinned by cross-multiplication (`x * x_den == x_num`),
  which absorbs the inverse application. The G2 tv2 inverses stay behind a
  single pair-inversion witness by default; `wide-witness` ships them
  separately (2 mul2 instead of 4, +96 B).
- Square-root witnesses travel in Montgomery form: the square check uses
  them as-is and reading the sign costs a bare reduction instead of the
  `to_mont` multiply.
- The min-pk verify fixture checks the pairing equation as one two-pair
  product, `e(pk, H) * e(-g1_gen, sig) == 1`: the second pair shares the
  final exponentiation a separate call pays again (one pair 25.7k CU, two
  pairs 38.7k), the negated generator is baked into the program instead of
  riding in the payload (96 fewer transaction bytes), and the GT identity
  compare uses the memcmp syscall (10 CU) instead of an inline 576 byte
  loop.
- The compact layouts (2026-07-13, tags 48/49/54/55): half-width square
  roots recovered through the curve equation, a single batched inverse for
  the whole pipeline, homogenized iso-3 denominators so the batch closes
  before any inverse exists, per-map iso with the E-side addition through
  the g2 add syscall, and an in-program binary-xgcd inverse replacing the
  batch witness in the 97 B layout. The real-part check runs as a
  difference of squares (one mont_mul against two mont_sqr) and the sign
  read is c0's parity (one from_mont against a full from_mont2). Extracting
  the shared iso-3 numerator and gx helpers moved the default G2 path
  253,404 to 252,916 (-489, better inlining); tag 51 followed to 308,209 at
  k=14.

## Measured dead ends

- 64-bit limb CIOS via u128: SBPF v3 dropped the v2-only PQR instruction
  class (UHMUL/SHMUL), so `(a as u128) * (b as u128)` lowers to a `__multi3`
  call at ~70 CU per multiply-accumulate (probe tag 21).
- Lazy Fp2 reduction (Karatsuba with unreduced 27-lane products and two wide
  reductions instead of three full multiplies): mul2 got ~9% slower. The
  fused multiply keeps operand splits shared and the accumulator hot;
  splitting product from reduction pays more in array traffic than the saved
  reduction pass (details at `fp2.rs::mul2`).
- Dedicated squaring under 32-bit CIOS did not pay; under the lane forms
  (ps29, now ps30) it does.
- An earlier SvdW direct-map variant (no isogeny) was dropped: it matches no
  standard suite and the witnessed SSWU pipeline now beats or ties it.

## Further optimizations

Open knobs, in rough order of interest:

- The modexp path (tags 30 to 33) runs hash_to_G1 in ~168k CU with zero
  witness bytes, against ~129k plus 338 B for the witnessed path. It needs
  `big_mod_exp` (SIMD-0529, merged but not active). Once 0529 activates, a
  transaction that is byte-bound rather than CU-bound should prefer it, and
  the compact G2 shape goes zero-witness too: the host `sqrt2` construction
  (an Fp2 root from three Fp modexps plus an inverse, the norm trick) runs
  through the syscall at ~1.7k CU per call, the branch decisions become
  Legendre modexps, and the ~112k xgcd inversion collapses to one modexp.
- The min-pk verify transaction is byte-bound, not CU-bound (308k to 638k
  of the 1.4M CU ceiling, but witness plus keys eat real transaction
  space). The byte/CU frontier is now the layout choice: 578 B at ~253k,
  145 B at ~470k, 97 B at ~583k, with `wide-witness` (674 B, ~241k) as the
  CU end once SIMD-0296 (4 KiB transactions, SDK support already merged)
  dissolves the byte constraint.
- `inv_xgcd` is deliberately naive: u and v stay full six-limb for the
  whole run, trailing-zero runs strip one bit per pass, and the odd-case
  halving materializes its add before shifting. Active-limb tracking plus
  a trailing-zeros batch is worth maybe 25-35k of its ~112k, at the cost
  of length-parameterized limb helpers in consensus-adjacent code.
- G2 cofactor clearing costs ~55k CU: roughly ~45k across ~140 g2 add syscalls
  and ~10k of psi/psi2 Fp2 multiplication. The Budroni-Pintore chain is the best
  known construction, so the syscall share is pricing rather than structure, but
  the endomorphism field work is real. The verify path feeds the hash into the
  pairing uncompressed, so no decompression cost hides there.
- The multiplier itself now runs ~1.93k CU for 338 products plus splits and
  packs. The remaining gap to the ALU floor is register pressure: 52 lanes
  live against 10 SBPF registers makes the allocator spill, and mul64
  destroying its destination charges a mov per cached operand. Fewer lanes
  would need a wider multiply, which the ISA does not have.
- The final validate syscall is defense-in-depth, not load-bearing: the
  cleared bytes come out of the last group-op syscall (on the curve by
  construction) and Budroni-Pintore style clearing lands any curve point in
  the subgroup. The pairing syscall subgroup-checks its own inputs (pinned
  by the syscall-contract test), so a pairing-bound consumer could drop the
  validate for another ~2.2k (G2) / ~1.8k (G1) CU; the standalone hash
  keeps it as the one runtime assertion on its own output. The stage
  parameter threaded through the prefix entry points costs ~3 CU end to
  end, measured tag against tag.
