//! BLS12-381 base field (Fp) arithmetic: product scanning Montgomery over
//! 30 bit lanes (ps30, R = 2^390 mod p) plus the big_mod_exp helpers used
//! by the host and the modexp assisted path.

use solana_program_error::ProgramError;

use crate::consts_g1::{INV, MODULUS, R2};
use crate::macros::{dot, lane, quotient};

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
    (s, c1 as u64 + c2 as u64)
}

#[inline(always)]
fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let (d, b1) = a.overflowing_sub(b);
    let (d, b2) = d.overflowing_sub(borrow);
    (d, b1 as u64 + b2 as u64)
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
/// digest any value below 2^384 (sums stay below 2p, products below
/// 2^768 < p 2^390, and the result is still below 2p)
#[inline(always)]
pub(crate) fn add_unreduced(a: &Fp, b: &Fp) -> Fp {
    let mut r = [0u64; 6];
    let mut carry = 0u64;
    for i in 0..6 {
        let (s, c) = adc(a[i], b[i], carry);
        r[i] = s;
        carry = c;
    }
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

/// Inverse of a Montgomery-form element, returned in Montgomery form.
pub(crate) fn inverse_mont(a: &Fp, exp_inv: &[u8; 48]) -> Result<Fp, ProgramError> {
    Ok(to_mont(&modexp(&from_mont(a), exp_inv)?))
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
    sum += spill;
    // column 11
    dot!(sum, a 0 1 2 3 4 5 6 7 8 9 10 11, b 11 10 9 8 7 6 5 4 3 2 1 0);
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10, P30 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 11);
    sum += spill;
    // column 12
    dot!(sum, a 0 1 2 3 4 5 6 7 8 9 10 11 12, b 12 11 10 9 8 7 6 5 4 3 2 1 0);
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10 11, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 12);
    sum += spill;
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
    sum += a[0] * a[0];
    quotient!(sum, m 0);
    // column 1
    dot!(sum, twice 0, a 1);
    dot!(sum, m 0, P30 1);
    quotient!(sum, m 1);
    // column 2
    dot!(sum, twice 0, a 2);
    sum += a[1] * a[1];
    dot!(sum, m 0 1, P30 2 1);
    quotient!(sum, m 2);
    // column 3
    dot!(sum, twice 0 1, a 3 2);
    dot!(sum, m 0 1 2, P30 3 2 1);
    quotient!(sum, m 3);
    // column 4
    dot!(sum, twice 0 1, a 4 3);
    sum += a[2] * a[2];
    dot!(sum, m 0 1 2 3, P30 4 3 2 1);
    quotient!(sum, m 4);
    // column 5
    dot!(sum, twice 0 1 2, a 5 4 3);
    dot!(sum, m 0 1 2 3 4, P30 5 4 3 2 1);
    quotient!(sum, m 5);
    // column 6
    dot!(sum, twice 0 1 2, a 6 5 4);
    sum += a[3] * a[3];
    dot!(sum, m 0 1 2 3 4 5, P30 6 5 4 3 2 1);
    quotient!(sum, m 6);
    // column 7
    dot!(sum, twice 0 1 2 3, a 7 6 5 4);
    dot!(sum, m 0 1 2 3 4 5 6, P30 7 6 5 4 3 2 1);
    quotient!(sum, m 7);
    // column 8
    dot!(sum, twice 0 1 2 3, a 8 7 6 5);
    sum += a[4] * a[4];
    dot!(sum, m 0 1 2 3 4 5 6 7, P30 8 7 6 5 4 3 2 1);
    quotient!(sum, m 8);
    // column 9
    dot!(sum, twice 0 1 2 3 4, a 9 8 7 6 5);
    dot!(sum, m 0 1 2 3 4 5 6 7 8, P30 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 9);
    // column 10
    dot!(sum, twice 0 1 2 3 4, a 10 9 8 7 6);
    sum += a[5] * a[5];
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9, P30 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 10);
    sum += spill;
    // column 11
    dot!(sum, twice 0 1 2 3 4 5, a 11 10 9 8 7 6);
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10, P30 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 11);
    sum += spill;
    // column 12
    dot!(sum, twice 0 1 2 3 4 5, a 12 11 10 9 8 7);
    sum += a[6] * a[6];
    let spill = sum >> 30;
    sum &= MASK30;
    dot!(sum, m 0 1 2 3 4 5 6 7 8 9 10 11, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    quotient!(sum, m 12);
    sum += spill;
    // column 13
    dot!(sum, twice 1 2 3 4 5 6, a 12 11 10 9 8 7);
    dot!(sum, m 1 2 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3 2 1);
    lane!(sum, r 0);
    // column 14
    dot!(sum, twice 2 3 4 5 6, a 12 11 10 9 8);
    sum += a[7] * a[7];
    dot!(sum, m 2 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3 2);
    lane!(sum, r 1);
    // column 15
    dot!(sum, twice 3 4 5 6 7, a 12 11 10 9 8);
    dot!(sum, m 3 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4 3);
    lane!(sum, r 2);
    // column 16
    dot!(sum, twice 4 5 6 7, a 12 11 10 9);
    sum += a[8] * a[8];
    dot!(sum, m 4 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5 4);
    lane!(sum, r 3);
    // column 17
    dot!(sum, twice 5 6 7 8, a 12 11 10 9);
    dot!(sum, m 5 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6 5);
    lane!(sum, r 4);
    // column 18
    dot!(sum, twice 6 7 8, a 12 11 10);
    sum += a[9] * a[9];
    dot!(sum, m 6 7 8 9 10 11 12, P30 12 11 10 9 8 7 6);
    lane!(sum, r 5);
    // column 19
    dot!(sum, twice 7 8 9, a 12 11 10);
    dot!(sum, m 7 8 9 10 11 12, P30 12 11 10 9 8 7);
    lane!(sum, r 6);
    // column 20
    dot!(sum, twice 8 9, a 12 11);
    sum += a[10] * a[10];
    dot!(sum, m 8 9 10 11 12, P30 12 11 10 9 8);
    lane!(sum, r 7);
    // column 21
    dot!(sum, twice 9 10, a 12 11);
    dot!(sum, m 9 10 11 12, P30 12 11 10 9);
    lane!(sum, r 8);
    // column 22
    dot!(sum, twice 10, a 12);
    sum += a[11] * a[11];
    dot!(sum, m 10 11 12, P30 12 11 10);
    lane!(sum, r 9);
    // column 23
    dot!(sum, twice 11, a 12);
    dot!(sum, m 11 12, P30 12 11);
    lane!(sum, r 10);
    // column 24
    sum += a[12] * a[12];
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
    sum += t[0];
    quotient!(sum, m 0);
    // column 1
    sum += t[1];
    dot!(sum, m 0, P30 1);
    quotient!(sum, m 1);
    // column 2
    sum += t[2];
    dot!(sum, m 0 1, P30 2 1);
    quotient!(sum, m 2);
    // column 3
    sum += t[3];
    dot!(sum, m 0 1 2, P30 3 2 1);
    quotient!(sum, m 3);
    // column 4
    sum += t[4];
    dot!(sum, m 0 1 2 3, P30 4 3 2 1);
    quotient!(sum, m 4);
    // column 5
    sum += t[5];
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
