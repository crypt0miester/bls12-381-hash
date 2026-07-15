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
| hash_to_G2 (RO, xgcd) | ~509k | 97 B | same suite, same output bytes |
| hash_to_G2 (RO, parity) | ~509k | 96 B | same suite, same output bytes |
| hash_to_G2 (RO, `wide-witness`) | ~241k | 674 B | same suite, bigger blob |
| encode_to_G1 (NU) | ~99k | 193 B | `_SSWU_NU_POP_`, byte-equal to blst encode |
| encode_to_G2 (NU) | ~173k | 385 B | `_SSWU_NU_POP_`, byte-equal to blst encode |

The default and `wide-witness` G2 rows are the same pipeline with a
different tv2 witness layout: the default 578 B blob pins both tv2 inverses
behind one pair inversion witness, while `wide-witness` ships both inverses
directly (96 more bytes, ~12k CU less) for when SIMD-0296 4 KiB
transactions land. The compact, xgcd and parity rows are the byte-bound
layouts, detailed below. An end-to-end min-pk BLS verify (hash_to_G2 plus
the pairing syscall) lands around 308k CU (297k wide, 525k compact, 564k
xgcd or parity). A 381-bit field multiplication costs ~1.9k CU with the
ps30 product-scanning multiplier (the textbook 32-bit CIOS form bottoms
out around 3.3k).

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

Three G2 layouts restructure the same checks to shrink the blob for
byte-bound transactions (hash tags 48/49/50, end-to-end verify tags
54/55/56, exact numbers at bench commit of 2026-07-14):

| layout | blob | hash CU | e2e verify k=14 / k=20 |
|---|---|---|---|
| default | 578 B | 252,914 | 308,220 / 294,567 |
| compact | 145 B | 470,152 | 525,426 / 511,815 |
| xgcd | 97 B | 509,080 | 564,353 / 550,742 |
| parity | 96 B | 509,029 | 564,303 / 550,692 |

The witness-free-inverse rows move a little with the message (the divsteps
batch count is input-dependent): over eight messages the parity hash spans
502.8k to 517.1k, average 512.0k. The 2024 hull bound caps any input at
1078 divsteps (36 batches before the g == 0 exit), putting the worst case
near ~519k.

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
batch product is inverted in-program by extended gcd
(`fp.rs::inv_divsteps`, Bernstein-Yang divsteps over 13 signed 30-bit
lanes; ~39k CU net of the witness check it replaces). Each batch of 30
divsteps runs on the low lanes and lands on the full values as two
row passes of signed multiply-adds; inside a batch, one multiple of f
cancels up to six low bits of g at a time through the odd-f identity
`f*(f^2-2) == -f^-1 mod 2^6`. The original theorem gives
floor((49*381 + 57)/17) = 1101 divsteps for 381-bit inputs (741 at 256
bits, the known number), which the 37-batch cap covers; the authors' 2024
convex-hull analysis tightens that to 1078, so at most 36 batches run
before the g == 0 exit (typically batch 27-28), and the tail requires
g == 0 and f == +-1 outright. `tools/check_divsteps.py` mirrors the code lane for
lane against `pow(a, p-2, p)`, proves the fused rounds equal the paper's
divstep, and tracks accumulator magnitudes (61 bits against the 63-bit
i64 budget); the host selftest re-checks products and the Fermat inverse
under debug overflow checks.

The parity layout then deletes the flags byte: the verifier
re-canonicalizes the root sign through sgn0 anyway, so which of the two
roots ships is a free bit, and each branch flag rides its root half's
parity (bit 0 of the last big-endian byte; `wit48` rejects non-canonical
encodings, so the bit is well defined, and p odd means the two roots
always differ in it). 96 B is the floor of this witness family: the two
c0 halves are pure computational advice for the square roots, and sqrt
mod p stays exponentiation-shaped with no multiply-free algorithm known.

A branch lie either desyncs the batch from `w` or leaves the real-part
check unsatisfiable (the wrong branch's gx is a non-square); the
soundness sweeps cover all layouts bit by bit, including
batch-consistent lies (`generate_compact_*_steered`) and the other
square root (`flip_first_root`, which also negates the trailing inverse
when the layout carries one). In the parity layout the other root IS a
branch lie, so the sweep asserts it aborts: each message has exactly one
acceptable blob, deleting witness malleability.

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
  the g2 add syscall, and an in-program extended-gcd inverse replacing the
  batch witness in the 97 B layout. The real-part check runs as a
  difference of squares (one mont_mul against two mont_sqr) and the sign
  read is c0's parity (one from_mont against a full from_mont2). Extracting
  the shared iso-3 numerator and gx helpers moved the default G2 path
  253,404 to 252,916 (-489, better inlining); tag 51 followed to 308,209 at
  k=14.
- Divsteps replaced the binary xgcd (2026-07-14): the bit-at-a-time
  shift-and-subtract inverse (~112k CU) became Bernstein-Yang divsteps in
  30-bit batches (~39k net). The winning shape on this ISA: branchy
  divsteps (a masked constant-time round measured 63 instructions per
  divstep against ~5 fused), w-rounds cancelling up to 6 bits of g per
  multiple of f, straight-line row passes (one matrix row live per pass,
  under the ten-register file), update calls kept `inline(never)` with
  ping-pong output buffers swapped by reference. An R3 constant
  (`consts_g1.rs`) returns the inverse to the Montgomery domain in one
  mont_mul, replacing the from_mont/to_mont pair around the old xgcd.
  Tags 49/55 moved 582,583 to 509,080 and 637,851 to 564,353 (k=14).
- The parity layout (2026-07-14, tags 50/56): the SSWU branch bits ride
  the root halves' parity instead of a flags byte, 97 B to 96 B for ~51
  CU, and the witness becomes unique per message (the other root reads as
  a branch lie and aborts).

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
- Constant-time divsteps: the fully masked branchless round (libsecp256k1's
  modinv32 shape) measured 63 instructions per divstep on sBPF, where a
  taken branch costs one instruction and the mask emulation of one
  conditional costs six, with the 13-value live set spilling on top.
  Branchy per-step form: ~18. Fused w-rounds: ~5. There is no secret data
  on-chain, so constant time buys nothing here.
- Inlining the divsteps update passes (or unrolling the batch loop by two
  to hard-code the ping-pong buffers): 39.5k against 34.8k for the
  `inline(never)` calls; the merged live sets spill more than the call
  overhead saves.
- Half-delta divsteps (delta starting at 1/2, the 2024 follow-up's
  variant): the proven bound drops 1078 to 878 divsteps at 381 bits, but
  the typical convergence barely moves (~2.1 divsteps per bit either
  way; model worst 840 to 810) while delta hugging zero shrinks the
  average w-round limit, so the fused rounds amortize worse: probe 36.1k
  against 34.8k, parity hash 513,371 against 509,029. The tighter cap
  only helps the worst-case bound, which the delta = 1 hull number (1078,
  36 batches) already serves. A per-step or masked implementation would
  benefit instead; this fused-vartime shape does not.
- 60-step composed divsteps batches (the 2024 follow-up's wide-update
  idea): two 30-step batches composed into 60-bit matrix entries and
  applied as one row pass each, with balanced 30-bit entry splits and a
  two-limb quotient. Verified correct (record and bounds in
  `tools/experiment_divsteps60.py`) and measured a wash: 34,590.6
  against 34,740.0 per inverse. Composition conserves the product count,
  so the win had to come from halved pass overhead, and the fused row
  pass carries four split entries, two quotient limbs, four sliding
  operands and an accumulator past the ten registers; the spills give
  the overhead back.
- Packing (f, u, v) and (g, q, r) into one word each for the inner loop
  (same talk): rejected on arithmetic before code. Fixed-width fields
  cannot hold a w-round (|w u| reaches 2^36), so the full sliding-
  boundary scheme is forced, and its unpack-compose overhead (~45 per
  batch) cancels most of the per-step savings (strips 8 to 5, w-rounds
  ~25 to ~14, swaps ~10 to 4): the optimistic ceiling is ~1.4k on an
  inner loop that is ~8k of a 509k pipeline, below the noise the two
  measured washes above landed in.

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
- The min-pk verify transaction is byte-bound, not CU-bound (308k to 564k
  of the 1.4M CU ceiling, but witness plus keys eat real transaction
  space). The byte/CU frontier is now the layout choice: 578 B at ~253k,
  145 B at ~470k, 96 B at ~509k, with `wide-witness` (674 B, ~241k) as the
  CU end once SIMD-0296 (4 KiB transactions, SDK support already merged)
  dissolves the byte constraint.
- `inv_divsteps` has maybe 4-6k of allocator slop left (the update passes
  execute ~30% more instructions than their arithmetic core; LLVM keeps
  spilling row entries), and dynamic f/g length tracking is worth ~2k
  more. Both were measured as small against the structural moves that
  landed, and the next factor would need hand-written sBPF, which this
  crate does not want. The batch geometry itself is optimal: total update
  work scales as 1/(batch_bits * lane_bits) under batch_bits + lane_bits
  <= 61, maximized at the current 30/30.
- The safegcd authors' December 2024 follow-up (Bernstein, Chen,
  Harrison, Maxwell, Wang, Wuille, Yang,
  [More on fast constant-time gcd computation and modular inversion](https://troll.iis.sinica.edu.tw/ws2024/safegcd2.pdf),
  with hull computations in
  [sipa/safegcd-bounds](https://github.com/sipa/safegcd-bounds)) supplies
  the 1078-divstep bound the worst-case number rests on. Its three
  transplantable ideas are all resolved: the half-delta variant and the
  60-step composed updates were measured and rejected, the packed inner
  loop rejected on arithmetic (all three in the dead ends).
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
