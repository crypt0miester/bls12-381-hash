//! BLS12-381 base field (Fp) arithmetic: product scanning Montgomery over
//! 30 bit lanes (ps30, R = 2^390 mod p) plus the big_mod_exp helpers used
//! by the host and the modexp assisted path.

use solana_program_error::ProgramError;

use crate::consts_g1::{INV, MODULUS, R2, R3};
use crate::macros::{dot, lane, quotient, row_limbs, row_limbs_fold, row_limbs_modp};

pub(crate) type Fp = [u64; 6];

// Split a tag payload into (witness, message), erroring if it is too short.
pub(crate) fn split_witness(payload: &[u8], total: usize) -> Result<(&[u8], &[u8]), ProgramError> {
    if payload.len() < total {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(payload.split_at(total))
}

const ONE: Fp = [1, 0, 0, 0, 0, 0];

#[inline(always)]
fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let (s, c1) = a.overflowing_add(b);
    let (s, c2) = s.overflowing_add(carry);
    (s, (c1 as u64).wrapping_add(c2 as u64))
}

#[inline(always)]
fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let (d, b1) = a.overflowing_sub(b);
    let (d, b2) = d.overflowing_sub(borrow);
    (d, (b1 as u64).wrapping_add(b2 as u64))
}

#[inline(always)]
pub(crate) fn geq(a: &Fp, b: &Fp) -> bool {
    for i in (0..6).rev() {
        if a[i] > b[i] {
            return true;
        }
        if a[i] < b[i] {
            return false;
        }
    }
    true
}

#[inline(always)]
pub(crate) fn sub_nocheck(a: &Fp, b: &Fp) -> Fp {
    let mut r = [0u64; 6];
    let mut borrow = 0u64;
    for i in 0..6 {
        let (d, br) = sbb(a[i], b[i], borrow);
        r[i] = d;
        borrow = br;
    }
    r
}

#[inline(always)]
pub(crate) fn add_mod(a: &Fp, b: &Fp) -> Fp {
    let mut r = [0u64; 6];
    let mut carry = 0u64;
    for i in 0..6 {
        let (s, c) = adc(a[i], b[i], carry);
        r[i] = s;
        carry = c;
    }
    if carry != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

#[inline(always)]
pub(crate) fn sub_mod(a: &Fp, b: &Fp) -> Fp {
    if geq(a, b) {
        sub_nocheck(a, b)
    } else {
        add_carryless(&sub_nocheck(a, b))
    }
}

// wrapped subtraction result plus p restores the field value
#[inline(always)]
pub(crate) fn add_carryless(r: &Fp) -> Fp {
    let mut out = [0u64; 6];
    let mut carry = 0u64;
    for i in 0..6 {
        let (s, c) = adc(r[i], MODULUS[i], carry);
        out[i] = s;
        carry = c;
    }
    out
}

#[inline(always)]
pub(crate) fn neg_mod(a: &Fp) -> Fp {
    if is_zero(a) {
        return [0u64; 6];
    }
    sub_nocheck(&MODULUS, a)
}

#[inline(always)]
pub(crate) fn is_zero(a: &Fp) -> bool {
    (a[0] | a[1] | a[2] | a[3] | a[4] | a[5]) == 0
}

/// a + b without the modular reduction; only for mont_mul operands, which
/// digest any value below 2^384.
#[inline(always)]
pub(crate) fn add_unreduced(a: &Fp, b: &Fp) -> Fp {
    let mut r = [0u64; 6];
    let mut carry = 0u64;
    for i in 0..6 {
        let (s, c) = adc(a[i], b[i], carry);
        r[i] = s;
        carry = c;
    }
    debug_assert_eq!(carry, 0, "add_unreduced: operand sum carried past 2^384");
    r
}

pub(crate) fn to_mont(a: &Fp) -> Fp {
    mont_mul(a, &R2)
}

pub(crate) fn limbs_to_be(a: &Fp) -> [u8; 48] {
    let mut out = [0u8; 48];
    for i in 0..6 {
        out[i * 8..i * 8 + 8].copy_from_slice(&a[5 - i].to_be_bytes());
    }
    out
}

pub(crate) fn be_to_limbs(b: &[u8; 48]) -> Fp {
    let mut r = [0u64; 6];
    for i in 0..6 {
        r[5 - i] = u64::from_be_bytes(b[i * 8..i * 8 + 8].try_into().unwrap());
    }
    r
}

pub(crate) fn shr1(a: &Fp) -> Fp {
    let mut r = [0u64; 6];
    for i in 0..6 {
        r[i] = a[i] >> 1;
        if i < 5 {
            r[i] |= a[i + 1] << 63;
        }
    }
    r
}

pub(crate) fn exp_inverse() -> [u8; 48] {
    limbs_to_be(&sub_nocheck(&MODULUS, &[2, 0, 0, 0, 0, 0]))
}

pub(crate) fn exp_legendre() -> [u8; 48] {
    limbs_to_be(&shr1(&sub_nocheck(&MODULUS, &ONE)))
}

pub(crate) fn exp_sqrt() -> [u8; 48] {
    let mut p1 = MODULUS;
    p1[0] += 1;
    limbs_to_be(&shr1(&shr1(&p1)))
}

#[cfg(target_os = "solana")]
pub(crate) mod sys {
    use solana_define_syscall::define_syscall;

    define_syscall!(fn sol_curve_validate_point(curve_id: u64, point_addr: *const u8, result: *mut u8) -> u64);
    define_syscall!(fn sol_curve_group_op(curve_id: u64, group_op: u64, left_input_addr: *const u8, right_input_addr: *const u8, result_point_addr: *mut u8) -> u64);
    define_syscall!(fn sol_big_mod_exp(params: *const u8, result: *mut u8) -> u64);
}

#[cfg(not(target_os = "solana"))]
#[allow(clippy::missing_safety_doc)]
pub(crate) mod sys {
    pub unsafe fn sol_curve_validate_point(_: u64, _: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_curve_group_op(_: u64, _: u64, _: *const u8, _: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_big_mod_exp(_: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
}

#[repr(C)]
struct BigModExpParams {
    base: *const u8,
    base_len: u64,
    exponent: *const u8,
    exponent_len: u64,
    modulus: *const u8,
    modulus_len: u64,
}

pub(crate) fn modexp_bytes(base: &[u8], exp: &[u8]) -> Result<[u8; 48], ProgramError> {
    let modulus = limbs_to_be(&MODULUS);
    let params = BigModExpParams {
        base: base.as_ptr(),
        base_len: base.len() as u64,
        exponent: exp.as_ptr(),
        exponent_len: exp.len() as u64,
        modulus: modulus.as_ptr(),
        modulus_len: 48,
    };
    // The syscall fills all 48 bytes on success, so skip the zero-init
    let mut out = core::mem::MaybeUninit::<[u8; 48]>::uninit();
    let rc = unsafe {
        sys::sol_big_mod_exp(
            &params as *const BigModExpParams as *const u8,
            out.as_mut_ptr() as *mut u8,
        )
    };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    // SAFETY: rc == 0 means the syscall wrote all 48 bytes
    Ok(unsafe { out.assume_init() })
}

/// Modular exponentiation of a canonical-form element, returning canonical.
pub(crate) fn modexp(base: &Fp, exp: &[u8; 48]) -> Result<Fp, ProgramError> {
    Ok(be_to_limbs(&modexp_bytes(&limbs_to_be(base), exp)?))
}

/// Inverse of a Montgomery-form element, returned in Montgomery form: the
/// Montgomery representative feeds the syscall as is ((aR)^(p-2) =
/// a^-1 R^-1) and the R3 multiply returns to the domain, skipping the
/// from_mont a canonical round trip would pay.
pub(crate) fn inverse_mont(a: &Fp, exp_inv: &[u8; 48]) -> Result<Fp, ProgramError> {
    Ok(mont_mul(&modexp(a, exp_inv)?, &R3))
}

#[inline(always)]
pub(crate) fn is_one(a: &Fp) -> bool {
    a[0] == 1 && (a[1] | a[2] | a[3] | a[4] | a[5]) == 0
}

/// Halve a residue mod p: even values shift, odd values add p first
/// (x1 < p keeps the sum under 2^382, so no carry is lost).
#[inline(always)]
pub(crate) fn half_mod(a: &Fp) -> Fp {
    if a[0] & 1 == 0 {
        shr1(a)
    } else {
        shr1(&add_carryless(a))
    }
}

// Modular inverse by Bernstein-Yang divsteps ("safegcd"): up to 37 batches
// of 30 divsteps on the low lanes, each batch applied to the full values
// as signed 13-lane row passes. The 2024 hull analysis bounds 381-bit
// inputs by 1078 divsteps (36 batches; the cap keeps the older 1101
// bound's margin) and the loop exits as soon as g == 0.
// tools/check_divsteps.py mirrors this code lane for lane and carries
// the bound and overflow arguments.

const M30S: i64 = MASK30 as i64;

/// p^-1 mod 2^30 (positive, cut from the negated 64-bit INV)
const PINV30: i64 = (INV.wrapping_neg() & MASK30) as i64;
const _: () = assert!((P30[0] as i64).wrapping_mul(PINV30) & M30S == 1);

/// The modulus in signed lanes for the d/e updates
const P30S: [i64; 13] = {
    let mut r = [0i64; 13];
    let mut i = 0;
    while i < 13 {
        r[i] = P30[i] as i64;
        i += 1;
    }
    r
};

/// Signed-lane values: lanes 0..11 proper in [0, 2^30), lane 12 carries
/// the sign; everything stays within (-2p, 2p), so lane 12 fits 23 bits.
type S13 = [i64; 13];

#[inline(always)]
fn split30_signed(x: &Fp) -> S13 {
    let u = split30(x);
    let mut r = [0i64; 13];
    let mut i = 0;
    while i < 13 {
        r[i] = u[i] as i64;
        i += 1;
    }
    r
}

/// 30 divsteps on the low lanes, batched, with eta = -delta. Branchy on
/// purpose: nothing on chain is secret, and a taken branch beats a masked
/// select on sBPF. Each w-round cancels min(eta+1, remaining, 6) low bits
/// of g with one multiple of f, using the odd-f identity
/// f*(f^2-2) == -f^-1 mod 2^6, and shifts them out in one go: the first
/// consumed step is the odd add, the rest are the even halvings, and an
/// odd add at consumed step j needs delta + j - 1 <= 0, hence the eta+1
/// cap. The f row doubles per consumed step instead of g halving, so the
/// matrix stays integral, t/2^30 is the true transition, and |u|+|v| and
/// |q|+|r| stay within 2^30. In-loop |f|, |g| stay under 2^32 and the
/// w-round transients under 2^37.
#[inline(always)]
fn divsteps30(mut eta: i64, f0: i64, g0: i64) -> (i64, i64, i64, i64, i64) {
    let (mut u, mut v, mut q, mut r) = (1i64, 0i64, 0i64, 1i64);
    let (mut f, mut g) = (f0, g0);
    let mut i = 30i64;
    loop {
        debug_assert!(f & 1 == 1);
        // strip trailing zeros of g one at a time, up to the step budget
        while g & 1 == 0 {
            g >>= 1;
            u <<= 1;
            v <<= 1;
            eta = eta.wrapping_sub(1);
            i -= 1;
            if i == 0 {
                return (eta, u, v, q, r);
            }
        }
        if eta < 0 {
            eta = eta.wrapping_neg();
            core::mem::swap(&mut f, &mut g);
            g = g.wrapping_neg();
            core::mem::swap(&mut u, &mut q);
            q = q.wrapping_neg();
            core::mem::swap(&mut v, &mut r);
            r = r.wrapping_neg();
        }
        let limit = eta.wrapping_add(1).min(i).min(6);
        let m = (1i64 << limit) - 1;
        let w = g.wrapping_mul(f.wrapping_mul(f.wrapping_mul(f).wrapping_sub(2))) & m;
        debug_assert!(w & 1 == 1);
        g = g.wrapping_add(w.wrapping_mul(f));
        debug_assert!(g & m == 0);
        q = q.wrapping_add(w.wrapping_mul(u));
        r = r.wrapping_add(w.wrapping_mul(v));
        g >>= limit;
        u <<= limit;
        v <<= limit;
        eta = eta.wrapping_sub(limit);
        i -= limit;
        if i == 0 {
            return (eta, u, v, q, r);
        }
    }
}

/// One batch applied to f and g: t*(f, g) / 2^30, the division exact by
/// the matrix construction. Two row passes so each keeps one matrix row
/// and one accumulator live: the f row lands in fout (the caller swaps
/// buffers), the g row updates in place, reading the untouched old f.
/// Returns whether g became zero (the or-fold is free next to the
/// stores), which ends the outer loop early: once g == 0 every further
/// batch is the identity on f, d and e (matrix [[2^30, 0], [0, 1]],
/// cancelled by the shared division).
#[inline(always)]
fn update_fg(f: &S13, g: &mut S13, fout: &mut S13, u: i64, v: i64, q: i64, r: i64) -> bool {
    let mut cf = u.wrapping_mul(f[0]).wrapping_add(v.wrapping_mul(g[0]));
    debug_assert!(cf & M30S == 0);
    cf >>= 30;
    row_limbs!(fout f g cf, u v; 1 2 3 4 5 6 7 8 9 10 11 12);
    fout[12] = cf;
    let mut cg = q.wrapping_mul(f[0]).wrapping_add(r.wrapping_mul(g[0]));
    debug_assert!(cg & M30S == 0);
    cg >>= 30;
    let mut nonzero = 0i64;
    row_limbs_fold!(g f g cg nonzero, q r; 1 2 3 4 5 6 7 8 9 10 11 12);
    g[12] = cg;
    (nonzero | cg) == 0
}

/// One batch applied to d and e mod p: t*(d, e) / 2^30 with the division
/// realized by adding the p-multiple that zeroes the low 30 bits (md/me,
/// one wrapping multiply by PINV30 each); the sign masks fold an extra p
/// in for negative d/e, keeping both in (-2p, p). Two row passes like
/// update_fg; the head quotients need the old d[0] and e[0], so both are
/// fixed before either pass stores.
#[inline(always)]
fn update_de(d: &S13, e: &mut S13, dout: &mut S13, u: i64, v: i64, q: i64, r: i64) {
    let sd = d[12] >> 63;
    let se = e[12] >> 63;
    let mut md = (u & sd).wrapping_add(v & se);
    let mut me = (q & sd).wrapping_add(r & se);
    let mut cd = u.wrapping_mul(d[0]).wrapping_add(v.wrapping_mul(e[0]));
    let mut ce = q.wrapping_mul(d[0]).wrapping_add(r.wrapping_mul(e[0]));
    md = md.wrapping_sub(PINV30.wrapping_mul(cd).wrapping_add(md) & M30S);
    me = me.wrapping_sub(PINV30.wrapping_mul(ce).wrapping_add(me) & M30S);
    cd = cd.wrapping_add(P30S[0].wrapping_mul(md));
    ce = ce.wrapping_add(P30S[0].wrapping_mul(me));
    debug_assert!(cd & M30S == 0 && ce & M30S == 0);
    cd >>= 30;
    ce >>= 30;
    row_limbs_modp!(dout d e cd, u v md; 1 2 3 4 5 6 7 8 9 10 11 12);
    dout[12] = cd;
    row_limbs_modp!(e d e ce, q r me; 1 2 3 4 5 6 7 8 9 10 11 12);
    e[12] = ce;
}

/// Conditionally negate d (in (-2p, p)), lift by p while negative (two
/// masked passes cover the worst case), and pack canonical
fn signed_to_fp(d: &S13, negate: i64) -> Fp {
    let mut t = [0i64; 13];
    let mut carry = 0i64;
    for i in 0..12 {
        let lane = ((d[i] ^ negate).wrapping_sub(negate)).wrapping_add(carry);
        t[i] = lane & M30S;
        carry = lane >> 30;
    }
    t[12] = ((d[12] ^ negate).wrapping_sub(negate)).wrapping_add(carry);
    for _ in 0..2 {
        let lift = t[12] >> 63;
        let mut carry = 0i64;
        for i in 0..12 {
            let lane = t[i].wrapping_add(P30S[i] & lift).wrapping_add(carry);
            t[i] = lane & M30S;
            carry = lane >> 30;
        }
        t[12] = t[12].wrapping_add((P30S[12] & lift).wrapping_add(carry));
    }
    debug_assert!(t[12] >= 0);
    let mut out = [0u64; 13];
    let mut i = 0;
    while i < 13 {
        out[i] = t[i] as u64;
        i += 1;
    }
    pack30(&out)
}

/// Inverse of a nonzero residue below p, of the value as given (feed it a
/// Montgomery-form element and multiply by R3 to return to that domain).
/// The cost moves mildly with the input through the batch count, bounded
/// by the cap.
pub(crate) fn inv_divsteps(a: &Fp) -> Result<Fp, ProgramError> {
    if is_zero(a) {
        return Err(ProgramError::InvalidInstructionData);
    }
    let mut d_a = [0i64; 13];
    let mut d_b = [0i64; 13];
    let mut e = [0i64; 13];
    e[0] = 1;
    let mut f_a = P30S;
    let mut f_b = [0i64; 13];
    let mut g = split30_signed(a);
    // the row passes write f and d into the spare buffer; swapping the
    // references costs a register move, not a copy
    let (mut f, mut f_next) = (&mut f_a, &mut f_b);
    let (mut d, mut d_next) = (&mut d_a, &mut d_b);
    let mut eta = -1i64;
    let mut done = false;
    for _ in 0..37 {
        let (eta_next, u, v, q, r) = divsteps30(eta, f[0], g[0]);
        eta = eta_next;
        done = update_fg(f, &mut g, f_next, u, v, q, r);
        update_de(d, &mut e, d_next, u, v, q, r);
        core::mem::swap(&mut f, &mut f_next);
        core::mem::swap(&mut d, &mut d_next);
        if done {
            break;
        }
    }
    // gcd(a, p) = 1 for 0 < a < p prime, so these are unreachable defense
    // in depth like the old xgcd pass cap
    let one: S13 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let minus_one: S13 = [
        M30S, M30S, M30S, M30S, M30S, M30S, M30S, M30S, M30S, M30S, M30S, M30S, -1,
    ];
    if !done || (*f != one && *f != minus_one) {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(signed_to_fp(d, f[12] >> 63))
}


pub(crate) fn wit48(bytes: &[u8]) -> Result<Fp, ProgramError> {
    let arr: &[u8; 48] = bytes
        .try_into()
        .map_err(|_| ProgramError::InvalidInstructionData)?;
    let limbs = be_to_limbs(arr);
    if geq(&limbs, &MODULUS) {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(limbs)
}


// Product scanning Montgomery multiplication over 30 bit lanes ("ps30").
//
// The bodies are generated by tools/gen_ps30.py and expand to straight
// line code. One lane fewer than the 29 bit form pays for the tall middle
// columns summing past 2^64: there the generated code banks the high bits
// into a spill register between the operand run and the modulus run. The
// generator tracks the exact worst case of every column sum and carry.

const MASK30: u64 = (1 << 30) - 1;

/// The negated modulus inverse mod 2^30, cut down from the 64 bit INV
const INV30: u64 = INV & MASK30;

/// The modulus in 30 bit lanes, little endian
const P30: [u64; 13] = split30(&MODULUS);

/// Split six limbs into thirteen 30 bit lanes, little endian; spelled out
/// so every shift stays a constant
#[inline(always)]
const fn split30(x: &Fp) -> [u64; 13] {
    [
        x[0] & MASK30,
        (x[0] >> 30) & MASK30,
        ((x[0] >> 60) | (x[1] << 4)) & MASK30,
        (x[1] >> 26) & MASK30,
        ((x[1] >> 56) | (x[2] << 8)) & MASK30,
        (x[2] >> 22) & MASK30,
        ((x[2] >> 52) | (x[3] << 12)) & MASK30,
        (x[3] >> 18) & MASK30,
        ((x[3] >> 48) | (x[4] << 16)) & MASK30,
        (x[4] >> 14) & MASK30,
        ((x[4] >> 44) | (x[5] << 20)) & MASK30,
        (x[5] >> 10) & MASK30,
        x[5] >> 40,
    ]
}

/// Pack thirteen result lanes back into six limbs and canonicalize
#[inline(always)]
fn pack30(r: &[u64; 13]) -> Fp {
    let mut out = [
        r[0] | (r[1] << 30) | (r[2] << 60),
        (r[2] >> 4) | (r[3] << 26) | (r[4] << 56),
        (r[4] >> 8) | (r[5] << 22) | (r[6] << 52),
        (r[6] >> 12) | (r[7] << 18) | (r[8] << 48),
        (r[8] >> 16) | (r[9] << 14) | (r[10] << 44),
        (r[10] >> 20) | (r[11] << 10) | (r[12] << 40),
    ];
    if geq(&out, &MODULUS) {
        out = sub_nocheck(&out, &MODULUS);
    }
    out
}

/// hi || lo == u + q p over the integers, for hi || lo a 64-byte
/// big-endian value, u < 2^390 and q < 2^136 as 17 big-endian bytes.
/// Column sums over 30-bit lanes: five q lanes against the thirteen
/// modulus lanes keep every column below 2^63. Straight-line like the
/// multipliers, from tools/gen_quotient.py.
pub(crate) fn quotient_matches(hi: &[u8; 32], lo: &[u8; 32], u: &Fp, q_be: &[u8]) -> bool {
    let mut limbs = [0u64; 8];
    for i in 0..4 {
        limbs[7 - i] = u64::from_be_bytes(hi[8 * i..8 * i + 8].try_into().unwrap());
        limbs[3 - i] = u64::from_be_bytes(lo[8 * i..8 * i + 8].try_into().unwrap());
    }
    let ql = [
        u64::from_be_bytes(q_be[9..17].try_into().unwrap()),
        u64::from_be_bytes(q_be[1..9].try_into().unwrap()),
        q_be[0] as u64,
    ];
    let x = [
        limbs[0] & MASK30,
        (limbs[0] >> 30) & MASK30,
        ((limbs[0] >> 60) | (limbs[1] << 4)) & MASK30,
        (limbs[1] >> 26) & MASK30,
        ((limbs[1] >> 56) | (limbs[2] << 8)) & MASK30,
        (limbs[2] >> 22) & MASK30,
        ((limbs[2] >> 52) | (limbs[3] << 12)) & MASK30,
        (limbs[3] >> 18) & MASK30,
        ((limbs[3] >> 48) | (limbs[4] << 16)) & MASK30,
        (limbs[4] >> 14) & MASK30,
        ((limbs[4] >> 44) | (limbs[5] << 20)) & MASK30,
        (limbs[5] >> 10) & MASK30,
        ((limbs[5] >> 40) | (limbs[6] << 24)) & MASK30,
        (limbs[6] >> 6) & MASK30,
        ((limbs[6] >> 36) | (limbs[7] << 28)) & MASK30,
        (limbs[7] >> 2) & MASK30,
        (limbs[7] >> 32) & MASK30,
        (limbs[7] >> 62) & MASK30,
    ];
    let q = [
        ql[0] & MASK30,
        (ql[0] >> 30) & MASK30,
        ((ql[0] >> 60) | (ql[1] << 4)) & MASK30,
        (ql[1] >> 26) & MASK30,
        ((ql[1] >> 56) | (ql[2] << 8)) & MASK30,
    ];
    let u = split30(u);
    let mut col = 0u64;
    let mut diff = 0u64;
    col = col.wrapping_add(u[0]);
    dot!(col, q 0, P30 0);
    diff |= (col & MASK30) ^ x[0];
    col >>= 30;
    col = col.wrapping_add(u[1]);
    dot!(col, q 0 1, P30 1 0);
    diff |= (col & MASK30) ^ x[1];
    col >>= 30;
    col = col.wrapping_add(u[2]);
    dot!(col, q 0 1 2, P30 2 1 0);
    diff |= (col & MASK30) ^ x[2];
    col >>= 30;
    col = col.wrapping_add(u[3]);
    dot!(col, q 0 1 2 3, P30 3 2 1 0);
    diff |= (col & MASK30) ^ x[3];
    col >>= 30;
    col = col.wrapping_add(u[4]);
    dot!(col, q 0 1 2 3 4, P30 4 3 2 1 0);
    diff |= (col & MASK30) ^ x[4];
    col >>= 30;
    col = col.wrapping_add(u[5]);
    dot!(col, q 0 1 2 3 4, P30 5 4 3 2 1);
    diff |= (col & MASK30) ^ x[5];
    col >>= 30;
    col = col.wrapping_add(u[6]);
    dot!(col, q 0 1 2 3 4, P30 6 5 4 3 2);
    diff |= (col & MASK30) ^ x[6];
    col >>= 30;
    col = col.wrapping_add(u[7]);
    dot!(col, q 0 1 2 3 4, P30 7 6 5 4 3);
    diff |= (col & MASK30) ^ x[7];
    col >>= 30;
    col = col.wrapping_add(u[8]);
    dot!(col, q 0 1 2 3 4, P30 8 7 6 5 4);
    diff |= (col & MASK30) ^ x[8];
    col >>= 30;
    col = col.wrapping_add(u[9]);
    dot!(col, q 0 1 2 3 4, P30 9 8 7 6 5);
    diff |= (col & MASK30) ^ x[9];
    col >>= 30;
    col = col.wrapping_add(u[10]);
    dot!(col, q 0 1 2 3 4, P30 10 9 8 7 6);
    diff |= (col & MASK30) ^ x[10];
    col >>= 30;
    col = col.wrapping_add(u[11]);
    dot!(col, q 0 1 2 3 4, P30 11 10 9 8 7);
    diff |= (col & MASK30) ^ x[11];
    col >>= 30;
    col = col.wrapping_add(u[12]);
    dot!(col, q 0 1 2 3 4, P30 12 11 10 9 8);
    diff |= (col & MASK30) ^ x[12];
    col >>= 30;
    dot!(col, q 1 2 3 4, P30 12 11 10 9);
    diff |= (col & MASK30) ^ x[13];
    col >>= 30;
    dot!(col, q 2 3 4, P30 12 11 10);
    diff |= (col & MASK30) ^ x[14];
    col >>= 30;
    dot!(col, q 3 4, P30 12 11);
    diff |= (col & MASK30) ^ x[15];
    col >>= 30;
    dot!(col, q 4, P30 12);
    diff |= (col & MASK30) ^ x[16];
    col >>= 30;
    diff |= (col & MASK30) ^ x[17];
    col >>= 30;
    (diff | col) == 0
}

/// x 2^-n mod p at compile time: odd values add p before the shift (x + p
/// stays below 2^382, so no carry leaves the top limb)
const fn c_halve(mut x: Fp, n: u32) -> Fp {
    let mut k = 0;
    while k < n {
        if x[0] & 1 == 1 {
            let mut carry = 0u64;
            let mut i = 0;
            while i < 6 {
                let (s, c1) = x[i].overflowing_add(MODULUS[i]);
                let (s, c2) = s.overflowing_add(carry);
                x[i] = s;
                carry = (c1 as u64) + (c2 as u64);
                i += 1;
            }
        }
        let mut i = 0;
        while i < 6 {
            x[i] = (x[i] >> 1) | if i < 5 { x[i + 1] << 63 } else { 0 };
            i += 1;
        }
        k += 1;
    }
    x
}

/// Lane tables for a fixed multiplier K < p: T[i] = K 2^(30 i - 330) mod p.
/// sum_i a_i T[i] is a K 2^60 R^-1 mod p up to multiples of p, so two
/// Montgomery quotient lanes (2^60) finish mont_mul(a, K) with 195 products
/// against 338, and every T lane rides as an immediate.
pub(crate) const fn fixed_table(k: Fp) -> [[u64; 13]; 13] {
    let mut t = [[0u64; 13]; 13];
    t[12] = split30(&crate::consts_g1::msn(k, 30));
    let mut v = k;
    let mut i = 12;
    while i > 0 {
        i -= 1;
        t[i] = split30(&v);
        v = c_halve(v, 30);
    }
    t
}

// Straight-line body of a fixed-multiplier product (tools/gen_fixed.py):
// column j sums the thirteen a_i T[i][j] plus the two quotient lanes'
// modulus terms, at most fifteen products below 2^60 and a carry below
// 2^35, so no column leaves the u64. The result stays below 2p.
macro_rules! fixed_mul {
    ($t:expr, $a:ident) => {{
        let a = split30($a);
        let mut r = [0u64; 13];
        let mut col = 0u64;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][0]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][0]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][0]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][0]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][0]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][0]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][0]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][0]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][0]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][0]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][0]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][0]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][0]));
        let m0 = col.wrapping_mul(INV30) & MASK30;
        col = col.wrapping_add(m0.wrapping_mul(P30[0])) >> 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][1]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][1]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][1]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][1]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][1]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][1]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][1]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][1]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][1]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][1]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][1]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][1]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][1]));
        col = col.wrapping_add(m0.wrapping_mul(P30[1]));
        let m1 = col.wrapping_mul(INV30) & MASK30;
        col = col.wrapping_add(m1.wrapping_mul(P30[0])) >> 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][2]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][2]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][2]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][2]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][2]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][2]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][2]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][2]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][2]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][2]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][2]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][2]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][2]));
        col = col.wrapping_add(m0.wrapping_mul(P30[2])).wrapping_add(m1.wrapping_mul(P30[1]));
        r[0] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][3]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][3]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][3]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][3]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][3]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][3]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][3]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][3]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][3]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][3]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][3]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][3]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][3]));
        col = col.wrapping_add(m0.wrapping_mul(P30[3])).wrapping_add(m1.wrapping_mul(P30[2]));
        r[1] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][4]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][4]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][4]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][4]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][4]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][4]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][4]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][4]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][4]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][4]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][4]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][4]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][4]));
        col = col.wrapping_add(m0.wrapping_mul(P30[4])).wrapping_add(m1.wrapping_mul(P30[3]));
        r[2] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][5]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][5]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][5]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][5]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][5]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][5]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][5]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][5]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][5]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][5]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][5]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][5]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][5]));
        col = col.wrapping_add(m0.wrapping_mul(P30[5])).wrapping_add(m1.wrapping_mul(P30[4]));
        r[3] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][6]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][6]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][6]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][6]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][6]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][6]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][6]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][6]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][6]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][6]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][6]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][6]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][6]));
        col = col.wrapping_add(m0.wrapping_mul(P30[6])).wrapping_add(m1.wrapping_mul(P30[5]));
        r[4] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][7]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][7]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][7]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][7]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][7]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][7]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][7]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][7]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][7]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][7]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][7]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][7]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][7]));
        col = col.wrapping_add(m0.wrapping_mul(P30[7])).wrapping_add(m1.wrapping_mul(P30[6]));
        r[5] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][8]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][8]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][8]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][8]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][8]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][8]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][8]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][8]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][8]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][8]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][8]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][8]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][8]));
        col = col.wrapping_add(m0.wrapping_mul(P30[8])).wrapping_add(m1.wrapping_mul(P30[7]));
        r[6] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][9]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][9]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][9]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][9]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][9]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][9]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][9]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][9]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][9]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][9]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][9]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][9]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][9]));
        col = col.wrapping_add(m0.wrapping_mul(P30[9])).wrapping_add(m1.wrapping_mul(P30[8]));
        r[7] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][10]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][10]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][10]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][10]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][10]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][10]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][10]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][10]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][10]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][10]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][10]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][10]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][10]));
        col = col.wrapping_add(m0.wrapping_mul(P30[10])).wrapping_add(m1.wrapping_mul(P30[9]));
        r[8] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][11]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][11]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][11]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][11]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][11]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][11]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][11]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][11]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][11]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][11]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][11]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][11]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][11]));
        col = col.wrapping_add(m0.wrapping_mul(P30[11])).wrapping_add(m1.wrapping_mul(P30[10]));
        r[9] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(a[0].wrapping_mul($t[0][12]));
        col = col.wrapping_add(a[1].wrapping_mul($t[1][12]));
        col = col.wrapping_add(a[2].wrapping_mul($t[2][12]));
        col = col.wrapping_add(a[3].wrapping_mul($t[3][12]));
        col = col.wrapping_add(a[4].wrapping_mul($t[4][12]));
        col = col.wrapping_add(a[5].wrapping_mul($t[5][12]));
        col = col.wrapping_add(a[6].wrapping_mul($t[6][12]));
        col = col.wrapping_add(a[7].wrapping_mul($t[7][12]));
        col = col.wrapping_add(a[8].wrapping_mul($t[8][12]));
        col = col.wrapping_add(a[9].wrapping_mul($t[9][12]));
        col = col.wrapping_add(a[10].wrapping_mul($t[10][12]));
        col = col.wrapping_add(a[11].wrapping_mul($t[11][12]));
        col = col.wrapping_add(a[12].wrapping_mul($t[12][12]));
        col = col.wrapping_add(m0.wrapping_mul(P30[12])).wrapping_add(m1.wrapping_mul(P30[11]));
        r[10] = col & MASK30;
        col >>= 30;
        col = col.wrapping_add(m1.wrapping_mul(P30[12]));
        r[11] = col & MASK30;
        r[12] = col >> 30;
        pack30(&r)
    }};
}


/// One fixed multiplier per constant, out of line: each body is ~700
/// instructions of immediates.
macro_rules! fixed_mul_fn {
    ($name:ident, $table:ident) => {
        #[inline(never)]
        pub(crate) fn $name(a: &Fp) -> Fp {
            fixed_mul!($table, a)
        }
    };
}


/// Seven-limb accumulator for linear combinations with small integer
/// factors: terms add in unreduced and one reduction lands the sum.
pub(crate) type Wide = [u64; 7];

const MASK32: u64 = 0xffff_ffff;

/// acc += a k for a < 2^384 and k < 2^16, over 32-bit halves (the ISA has
/// no wide multiply)
#[inline(always)]
pub(crate) fn wide_mac(acc: &mut Wide, a: &Fp, k: u64) {
    let mut cm = 0u64;
    let mut ca = 0u64;
    for i in 0..6 {
        let lo = (a[i] & MASK32).wrapping_mul(k).wrapping_add(cm);
        let hi = (a[i] >> 32).wrapping_mul(k).wrapping_add(lo >> 32);
        cm = hi >> 32;
        let (s1, o1) = acc[i].overflowing_add((lo & MASK32) | (hi << 32));
        let (s2, o2) = s1.overflowing_add(ca);
        acc[i] = s2;
        ca = (o1 | o2) as u64;
    }
    acc[6] = acc[6].wrapping_add(cm).wrapping_add(ca);
}

/// acc += a
#[inline(always)]
pub(crate) fn wide_add(acc: &mut Wide, a: &Fp) {
    let mut c = 0u64;
    for i in 0..6 {
        let (s, c1) = adc(acc[i], a[i], c);
        acc[i] = s;
        c = c1;
    }
    acc[6] = acc[6].wrapping_add(c);
}

/// (p >> 336) + 1: dividing v >> 336 by it undershoots floor(v / p) by at
/// most one for any v below 2^400
const P_TOP_PLUS1: u64 = (MODULUS[5] >> 16) + 1;

/// v mod p for v < 2^400: one div64 estimates the quotient from the top
/// bits (never above floor(v / p), at most one below, since v >> 336 is
/// under 2^64 and p >> 336 over 2^44), a 32-bit-half multiply subtracts
/// it, and one conditional subtraction finishes.
#[inline(always)]
pub(crate) fn wide_reduce(v: &Wide) -> Fp {
    let q = ((v[5] >> 16) | (v[6] << 48)) / P_TOP_PLUS1;
    let mut r = [0u64; 6];
    let mut cm = 0u64;
    let mut br = 0u64;
    for i in 0..6 {
        let lo = (MODULUS[i] & MASK32).wrapping_mul(q).wrapping_add(cm);
        let hi = (MODULUS[i] >> 32).wrapping_mul(q).wrapping_add(lo >> 32);
        cm = hi >> 32;
        let (d, b) = sbb(v[i], (lo & MASK32) | (hi << 32), br);
        r[i] = d;
        br = b;
    }
    debug_assert_eq!(v[6].wrapping_sub(cm).wrapping_sub(br), 0);
    if geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    debug_assert!(!geq(&r, &MODULUS));
    r
}

/// Each operand lane doubled, for the squaring cross products
#[inline(always)]
fn double_lanes(a: &[u64; 13]) -> [u64; 13] {
    [
        a[0] << 1, a[1] << 1, a[2] << 1, a[3] << 1, a[4] << 1, a[5] << 1, a[6] << 1,
        a[7] << 1, a[8] << 1, a[9] << 1, a[10] << 1, a[11] << 1, a[12] << 1,
    ]
}

pub(crate) fn mont_mul(a: &Fp, b: &Fp) -> Fp {
    let a = split30(a);
    let b = split30(b);
    let mut m = [0u64; 13];
    let mut r = [0u64; 13];

    let mut sum = 0u64;
    // column 0
    dot!(sum, a 0, b 0);
    quotient!(sum, m 0);
    // column 1
    dot!(sum, a 0 1, b 1 0);
    dot!(sum, m 0, P30 1);
    quotient!(sum, m 1);
    // column 2
    dot!(sum, a 0 1 2, b 2 1 0);
    dot!(sum, m 0 1, P30 2 1);
    quotient!(sum, m 2);
    // column 3
    dot!(sum, a 0 1 2 3, b 3 2 1 0);
    dot!(sum, m 0 1 2, P30 3 2 1);
    quotient!(sum, m 3);
    // column 4
    dot!(sum, a 0 1 2 3 4, b 4 3 2 1 0);
    dot!(sum, m 0 1 2 3, P30 4 3 2 1);
    quotient!(sum, m 4);
    // column 5
    dot!(sum, a 0 1 2 3 4 5, b 5 4 3 2 1 0);
    dot!(sum, m 0 1 2 3 4, P30 5 4 3 2 1);
    quotient!(sum, m 5);
    // column 6
    dot!(sum, a 0 1 2 3 4 5 6, b 6 5 4 3 2 1 0);
    dot!(sum, m 0 1 2 3 4 5, P30 6 5 4 3 2 1);
    quotient!(sum, m 6);
    // column 7
    dot!(sum, a 0 1 2 3 4 5 6 7, b 7 6 5 4 3 2 1 0);
    dot!(sum, m 0 1 2 3 4 5 6, P30 7 6 5 4 3 2 1);
    quotient!(sum, m 7);
    // column 8
    dot!(sum, a 0 1 2 3 4 5 6 7 8, b 8 7 6 5 4 3 2 1 0);
    dot!(sum, m 0 1 2 3 4 5 6 7, P30 8 7 6 5 4 3 2 1);
    quotient!(sum, m 8);
    // column 9
    dot!(sum, a 0 1 2 3 4 5 6 7 8 9, b 9 8 7 6 5 4 3 2 1 0);
    dot!(sum, m 0 1 2 3 4 5 6 7 8, P30 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 9);
    // column 10
    dot!(sum, a 0 1 2 3 4 5 6 7 8 9 10, b 10 9 8 7 6 5 4 3 2 1 0);
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9, P30 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 10);
    sum = sum.wrapping_add(spill);
    // column 11
    dot!(sum, a 0 1 2 3 4 5 6 7 8 9 10 11, b 11 10 9 8 7 6 5 4 3 2 1 0);
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10, P30 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 11);
    sum = sum.wrapping_add(spill);
    // column 12
    dot!(sum, a 0 1 2 3 4 5 6 7 8 9 10 11 12, b 12 11 10 9 8 7 6 5 4 3 2 1 0);
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10 11, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 12);
    sum = sum.wrapping_add(spill);
    // column 13
    dot!(sum, a 1 2 3 4 5 6 7 8 9 10 11 12, b 12 11 10 9 8 7 6 5 4 3 2 1);
    dot!(sum, m 1 2 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    lane!(sum, r 0);
    // column 14
    dot!(sum, a 2 3 4 5 6 7 8 9 10 11 12, b 12 11 10 9 8 7 6 5 4 3 2);
    dot!(sum, m 2 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3 2);
    lane!(sum, r 1);
    // column 15
    dot!(sum, a 3 4 5 6 7 8 9 10 11 12, b 12 11 10 9 8 7 6 5 4 3);
    dot!(sum, m 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3);
    lane!(sum, r 2);
    // column 16
    dot!(sum, a 4 5 6 7 8 9 10 11 12, b 12 11 10 9 8 7 6 5 4);
    dot!(sum, m 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4);
    lane!(sum, r 3);
    // column 17
    dot!(sum, a 5 6 7 8 9 10 11 12, b 12 11 10 9 8 7 6 5);
    dot!(sum, m 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5);
    lane!(sum, r 4);
    // column 18
    dot!(sum, a 6 7 8 9 10 11 12, b 12 11 10 9 8 7 6);
    dot!(sum, m 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6);
    lane!(sum, r 5);
    // column 19
    dot!(sum, a 7 8 9 10 11 12, b 12 11 10 9 8 7);
    dot!(sum, m 7 8 9 10 11 12, P30 12 11 10 9 8 7);
    lane!(sum, r 6);
    // column 20
    dot!(sum, a 8 9 10 11 12, b 12 11 10 9 8);
    dot!(sum, m 8 9 10 11 12, P30 12 11 10 9 8);
    lane!(sum, r 7);
    // column 21
    dot!(sum, a 9 10 11 12, b 12 11 10 9);
    dot!(sum, m 9 10 11 12, P30 12 11 10 9);
    lane!(sum, r 8);
    // column 22
    dot!(sum, a 10 11 12, b 12 11 10);
    dot!(sum, m 10 11 12, P30 12 11 10);
    lane!(sum, r 9);
    // column 23
    dot!(sum, a 11 12, b 12 11);
    dot!(sum, m 11 12, P30 12 11);
    lane!(sum, r 10);
    // column 24
    dot!(sum, a 12, b 12);
    dot!(sum, m 12, P30 12);
    lane!(sum, r 11);
    debug_assert!(sum <= MASK30);
    r[12] = sum;
    pack30(&r)
}

pub(crate) fn mont_sqr(a: &Fp) -> Fp {
    let a = split30(a);
    let twice = double_lanes(&a);
    let mut m = [0u64; 13];
    let mut r = [0u64; 13];

    let mut sum = 0u64;
    // column 0
    sum = sum.wrapping_add(a[0].wrapping_mul(a[0]));
    quotient!(sum, m 0);
    // column 1
    dot!(sum, twice 0, a 1);
    dot!(sum, m 0, P30 1);
    quotient!(sum, m 1);
    // column 2
    dot!(sum, twice 0, a 2);
    sum = sum.wrapping_add(a[1].wrapping_mul(a[1]));
    dot!(sum, m 0 1, P30 2 1);
    quotient!(sum, m 2);
    // column 3
    dot!(sum, twice 0 1, a 3 2);
    dot!(sum, m 0 1 2, P30 3 2 1);
    quotient!(sum, m 3);
    // column 4
    dot!(sum, twice 0 1, a 4 3);
    sum = sum.wrapping_add(a[2].wrapping_mul(a[2]));
    dot!(sum, m 0 1 2 3, P30 4 3 2 1);
    quotient!(sum, m 4);
    // column 5
    dot!(sum, twice 0 1 2, a 5 4 3);
    dot!(sum, m 0 1 2 3 4, P30 5 4 3 2 1);
    quotient!(sum, m 5);
    // column 6
    dot!(sum, twice 0 1 2, a 6 5 4);
    sum = sum.wrapping_add(a[3].wrapping_mul(a[3]));
    dot!(sum, m 0 1 2 3 4 5, P30 6 5 4 3 2 1);
    quotient!(sum, m 6);
    // column 7
    dot!(sum, twice 0 1 2 3, a 7 6 5 4);
    dot!(sum, m 0 1 2 3 4 5 6, P30 7 6 5 4 3 2 1);
    quotient!(sum, m 7);
    // column 8
    dot!(sum, twice 0 1 2 3, a 8 7 6 5);
    sum = sum.wrapping_add(a[4].wrapping_mul(a[4]));
    dot!(sum, m 0 1 2 3 4 5 6 7, P30 8 7 6 5 4 3 2 1);
    quotient!(sum, m 8);
    // column 9
    dot!(sum, twice 0 1 2 3 4, a 9 8 7 6 5);
    dot!(sum, m 0 1 2 3 4 5 6 7 8, P30 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 9);
    // column 10
    dot!(sum, twice 0 1 2 3 4, a 10 9 8 7 6);
    sum = sum.wrapping_add(a[5].wrapping_mul(a[5]));
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9, P30 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 10);
    sum = sum.wrapping_add(spill);
    // column 11
    dot!(sum, twice 0 1 2 3 4 5, a 11 10 9 8 7 6);
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10, P30 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 11);
    sum = sum.wrapping_add(spill);
    // column 12
    dot!(sum, twice 0 1 2 3 4 5, a 12 11 10 9 8 7);
    sum = sum.wrapping_add(a[6].wrapping_mul(a[6]));
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10 11, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 12);
    sum = sum.wrapping_add(spill);
    // column 13
    dot!(sum, twice 1 2 3 4 5 6, a 12 11 10 9 8 7);
    dot!(sum, m 1 2 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    lane!(sum, r 0);
    // column 14
    dot!(sum, twice 2 3 4 5 6, a 12 11 10 9 8);
    sum = sum.wrapping_add(a[7].wrapping_mul(a[7]));
    dot!(sum, m 2 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3 2);
    lane!(sum, r 1);
    // column 15
    dot!(sum, twice 3 4 5 6 7, a 12 11 10 9 8);
    dot!(sum, m 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3);
    lane!(sum, r 2);
    // column 16
    dot!(sum, twice 4 5 6 7, a 12 11 10 9);
    sum = sum.wrapping_add(a[8].wrapping_mul(a[8]));
    dot!(sum, m 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4);
    lane!(sum, r 3);
    // column 17
    dot!(sum, twice 5 6 7 8, a 12 11 10 9);
    dot!(sum, m 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5);
    lane!(sum, r 4);
    // column 18
    dot!(sum, twice 6 7 8, a 12 11 10);
    sum = sum.wrapping_add(a[9].wrapping_mul(a[9]));
    dot!(sum, m 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6);
    lane!(sum, r 5);
    // column 19
    dot!(sum, twice 7 8 9, a 12 11 10);
    dot!(sum, m 7 8 9 10 11 12, P30 12 11 10 9 8 7);
    lane!(sum, r 6);
    // column 20
    dot!(sum, twice 8 9, a 12 11);
    sum = sum.wrapping_add(a[10].wrapping_mul(a[10]));
    dot!(sum, m 8 9 10 11 12, P30 12 11 10 9 8);
    lane!(sum, r 7);
    // column 21
    dot!(sum, twice 9 10, a 12 11);
    dot!(sum, m 9 10 11 12, P30 12 11 10 9);
    lane!(sum, r 8);
    // column 22
    dot!(sum, twice 10, a 12);
    sum = sum.wrapping_add(a[11].wrapping_mul(a[11]));
    dot!(sum, m 10 11 12, P30 12 11 10);
    lane!(sum, r 9);
    // column 23
    dot!(sum, twice 11, a 12);
    dot!(sum, m 11 12, P30 12 11);
    lane!(sum, r 10);
    // column 24
    sum = sum.wrapping_add(a[12].wrapping_mul(a[12]));
    dot!(sum, m 12, P30 12);
    lane!(sum, r 11);
    debug_assert!(sum <= MASK30);
    r[12] = sum;
    pack30(&r)
}

/// Out of Montgomery form: a times R^-1 mod p, the bare reduction
pub(crate) fn from_mont(x: &Fp) -> Fp {
    let t = split30(x);
    let mut m = [0u64; 13];
    let mut r = [0u64; 13];

    let mut sum = 0u64;
    // column 0
    sum = sum.wrapping_add(t[0]);
    quotient!(sum, m 0);
    // column 1
    sum = sum.wrapping_add(t[1]);
    dot!(sum, m 0, P30 1);
    quotient!(sum, m 1);
    // column 2
    sum = sum.wrapping_add(t[2]);
    dot!(sum, m 0 1, P30 2 1);
    quotient!(sum, m 2);
    // column 3
    sum = sum.wrapping_add(t[3]);
    dot!(sum, m 0 1 2, P30 3 2 1);
    quotient!(sum, m 3);
    // column 4
    sum = sum.wrapping_add(t[4]);
    dot!(sum, m 0 1 2 3, P30 4 3 2 1);
    quotient!(sum, m 4);
    // column 5
    sum = sum.wrapping_add(t[5]);
    dot!(sum, m 0 1 2 3 4, P30 5 4 3 2 1);
    quotient!(sum, m 5);
    // column 6
    sum += t[6];
    dot!(sum, m 0 1 2 3 4 5, P30 6 5 4 3 2 1);
    quotient!(sum, m 6);
    // column 7
    sum += t[7];
    dot!(sum, m 0 1 2 3 4 5 6, P30 7 6 5 4 3 2 1);
    quotient!(sum, m 7);
    // column 8
    sum += t[8];
    dot!(sum, m 0 1 2 3 4 5 6 7, P30 8 7 6 5 4 3 2 1);
    quotient!(sum, m 8);
    // column 9
    sum += t[9];
    dot!(sum, m 0 1 2 3 4 5 6 7 8, P30 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 9);
    // column 10
    sum += t[10];
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9, P30 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 10);
    // column 11
    sum += t[11];
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10, P30 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 11);
    // column 12
    sum += t[12];
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10 11, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 12);
    // column 13
    dot!(sum, m 1 2 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    lane!(sum, r 0);
    // column 14
    dot!(sum, m 2 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3 2);
    lane!(sum, r 1);
    // column 15
    dot!(sum, m 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3);
    lane!(sum, r 2);
    // column 16
    dot!(sum, m 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4);
    lane!(sum, r 3);
    // column 17
    dot!(sum, m 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5);
    lane!(sum, r 4);
    // column 18
    dot!(sum, m 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6);
    lane!(sum, r 5);
    // column 19
    dot!(sum, m 7 8 9 10 11 12, P30 12 11 10 9 8 7);
    lane!(sum, r 6);
    // column 20
    dot!(sum, m 8 9 10 11 12, P30 12 11 10 9 8);
    lane!(sum, r 7);
    // column 21
    dot!(sum, m 9 10 11 12, P30 12 11 10 9);
    lane!(sum, r 8);
    // column 22
    dot!(sum, m 10 11 12, P30 12 11 10);
    lane!(sum, r 9);
    // column 23
    dot!(sum, m 11 12, P30 12 11);
    lane!(sum, r 10);
    // column 24
    dot!(sum, m 12, P30 12);
    lane!(sum, r 11);
    debug_assert!(sum <= MASK30);
    r[12] = sum;
    pack30(&r)
}

/// Host guard for the small linear combinations: wide_mac / wide_add /
/// wide_reduce against add_mod chains, over edges and varied residues, up
/// to the largest combination the pipeline builds (566 p).
#[cfg(not(target_os = "solana"))]
pub fn wide_selftest() {
    use crate::consts_g1::{R, R2};
    let small = |a: &Fp, k: u64| {
        let mut r = [0u64; 6];
        for _ in 0..k {
            r = add_mod(&r, a);
        }
        r
    };
    let mut samples: alloc::vec::Vec<Fp> = alloc::vec![
        [0; 6],
        [1, 0, 0, 0, 0, 0],
        R,
        R2,
        sub_nocheck(&MODULUS, &[1, 0, 0, 0, 0, 0]),
    ];
    let mut x = R2;
    for i in 0..200u64 {
        x = add_mod(&mont_mul(&x, &R2), &[i, 3, 0, 0, 0, 0]);
        samples.push(x);
    }
    let top = sub_nocheck(&MODULUS, &[1, 0, 0, 0, 0, 0]);
    for a in &samples {
        for &(k1, k2, k3) in &[(1u64, 1u64, 1u64), (60, 253, 253), (6, 1, 1), (12, 1, 1), (1, 0, 0)] {
            let mut w = [0u64; 7];
            wide_mac(&mut w, a, k1);
            wide_mac(&mut w, &top, k2);
            wide_mac(&mut w, &neg_mod(a), k3);
            wide_add(&mut w, a);
            let want = add_mod(
                &add_mod(&small(a, k1), &small(&top, k2)),
                &add_mod(&small(&neg_mod(a), k3), a),
            );
            assert_eq!(wide_reduce(&w), want, "wide combination != add_mod chain");
        }
    }
}

/// The fixed multipliers the fat pipeline uses: each one equals
/// mont_mul(a, K) for its constant K, for any a below 2^384.
pub(crate) mod fixed {
    use super::*;
    use crate::consts_g2::{ISO3V_A3, ISO3V_CA3_HALF, PSI2_X_C0, PSI_X_C1, PSI_Y};

    const T_ISO3V_A3: [[u64; 13]; 13] = fixed_table(ISO3V_A3);
    const T_ISO3V_CA3_HALF: [[u64; 13]; 13] = fixed_table(ISO3V_CA3_HALF);
    const T_PSI_X_C1: [[u64; 13]; 13] = fixed_table(PSI_X_C1);
    const T_PSI_Y: [[u64; 13]; 13] = fixed_table(PSI_Y[0]);
    const T_PSI2_X_C0: [[u64; 13]; 13] = fixed_table(PSI2_X_C0);

    fixed_mul_fn!(iso3v_a3, T_ISO3V_A3);
    fixed_mul_fn!(iso3v_ca3_half, T_ISO3V_CA3_HALF);
    fixed_mul_fn!(psi_x_c1, T_PSI_X_C1);
    fixed_mul_fn!(psi_y, T_PSI_Y);
    fixed_mul_fn!(psi2_x_c0, T_PSI2_X_C0);

    /// Host guard: every fixed multiplier against mont_mul over edges,
    /// unreduced operands up to 2^384 and a chain of varied residues.
    #[cfg(not(target_os = "solana"))]
    pub fn selftest() {
        use crate::consts_g1::{R, R2};
        let cases: [(fn(&Fp) -> Fp, Fp); 5] = [
            (iso3v_a3, ISO3V_A3),
            (iso3v_ca3_half, ISO3V_CA3_HALF),
            (psi_x_c1, PSI_X_C1),
            (psi_y, PSI_Y[0]),
            (psi2_x_c0, PSI2_X_C0),
        ];
        let mut samples: alloc::vec::Vec<Fp> = alloc::vec![
            [0; 6],
            [1, 0, 0, 0, 0, 0],
            R,
            R2,
            sub_nocheck(&MODULUS, &[1, 0, 0, 0, 0, 0]),
            [u64::MAX; 6],
        ];
        let mut x = R2;
        for i in 0..300u64 {
            x = add_mod(&mont_mul(&x, &R2), &[i, 7, 0, 0, 0, 0]);
            samples.push(x);
            samples.push(add_unreduced(&add_unreduced(&x, &MODULUS), &MODULUS));
        }
        for (f, k) in cases {
            for s in &samples {
                assert_eq!(f(s), mont_mul(s, &k), "fixed multiplier != mont_mul");
            }
        }
    }
}
