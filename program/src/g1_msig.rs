//! Syscall-assisted RFC 9380 hash_to_G1 for BLS12-381 (min-sig suite
//! BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_).
//!
//! Strategy: expand_message_xmd through the sha256 syscall, wide reduction /
//! inversion / Legendre / sqrt through big_mod_exp, SSWU + iso-11 polynomial
//! evaluation with an in-program variable-time Montgomery multiplier, the two
//! mapped points added on the isogenous curve so the isogeny runs once, and
//! cofactor clearing as a double-and-add chain over the g1 add syscall.

use solana_program::program_error::ProgramError;

use crate::g1_consts::{
    INV, ISO11A_XDEN, ISO11A_XNUM, ISO11A_YDEN, ISO11A_YNUM, ISO11_XDEN, ISO11_XNUM, ISO11_YDEN,
    ISO11_YNUM, MODULUS, R, R2, SSWU_ELLP_A, SSWU_ELLP_B,
};
use crate::g2_consts::C256_MONT;

const DST_G1: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_";
const H_EFF: u64 = 0xd201000000010001;

pub(crate) type Fp = [u64; 6];

const ONE: Fp = [1, 0, 0, 0, 0, 0];

#[inline(always)]
fn mul_64x64(a: u64, b: u64) -> (u64, u64) {
    let a_lo = a & 0xffff_ffff;
    let a_hi = a >> 32;
    let b_lo = b & 0xffff_ffff;
    let b_hi = b >> 32;
    let p0 = a_lo * b_lo;
    let p1 = a_lo * b_hi;
    let p2 = a_hi * b_lo;
    let p3 = a_hi * b_hi;
    let mid = (p0 >> 32) + (p1 & 0xffff_ffff) + (p2 & 0xffff_ffff);
    let lo = (p0 & 0xffff_ffff) | (mid << 32);
    let hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
    (lo, hi)
}

#[inline(always)]
fn mac(acc: u64, a: u64, b: u64, carry: u64) -> (u64, u64) {
    let (lo, hi) = mul_64x64(a, b);
    let (lo, c1) = lo.overflowing_add(acc);
    let (lo, c2) = lo.overflowing_add(carry);
    (lo, hi + c1 as u64 + c2 as u64)
}

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

pub(crate) fn sub_mod(a: &Fp, b: &Fp) -> Fp {
    if geq(a, b) {
        sub_nocheck(a, b)
    } else {
        add_carryless(&sub_nocheck(a, b))
    }
}

pub(crate) fn add_carryless(r: &Fp) -> Fp {
    // wrapped subtraction result plus p restores the field value
    let mut out = [0u64; 6];
    let mut carry = 0u64;
    for i in 0..6 {
        let (s, c) = adc(r[i], MODULUS[i], carry);
        out[i] = s;
        carry = c;
    }
    out
}

pub(crate) fn neg_mod(a: &Fp) -> Fp {
    if a.iter().all(|&l| l == 0) {
        return [0u64; 6];
    }
    sub_nocheck(&MODULUS, a)
}

pub(crate) fn is_zero(a: &Fp) -> bool {
    a.iter().all(|&l| l == 0)
}

pub(crate) fn mont_mul(a: &Fp, b: &Fp) -> Fp {
    mont_mul_cios32(a, b)
}

/// Squaring entry point. the specialized SOS squaring 
/// below measures SLOWER on sBPF (4.3k vs 3.4k CU)
pub(crate) fn mont_sqr(a: &Fp) -> Fp {
    mont_mul_cios32(a, a)
}

/// Squaring with 32-bit limbs: cross products computed once and doubled,
/// then a standard Montgomery reduction pass. 
fn mont_sqr_sos32(a: &Fp) -> Fp {
    let mut a32 = [0u64; 12];
    let mut p32 = [0u64; 12];
    for i in 0..6 {
        a32[i * 2] = a[i] & 0xffff_ffff;
        a32[i * 2 + 1] = a[i] >> 32;
        p32[i * 2] = MODULUS[i] & 0xffff_ffff;
        p32[i * 2 + 1] = MODULUS[i] >> 32;
    }

    // 32-bit lanes of the 768-bit square, kept normalized per row
    let mut t = [0u64; 25];
    for i in 0..12 {
        let ai = a32[i];
        let mut carry = 0u64;
        for j in (i + 1)..12 {
            let (lo, hi) = mac32(t[i + j], ai, a32[j], carry);
            t[i + j] = lo;
            carry = hi;
        }
        let mut k = i + 12;
        let mut c = carry;
        while c != 0 && k < 24 {
            let v = t[k] + c;
            t[k] = v & 0xffff_ffff;
            c = v >> 32;
            k += 1;
        }
        t[24] += c;
    }

    // double the cross products
    let mut c = 0u64;
    for lane in t.iter_mut().take(24) {
        let v = (*lane << 1) | c;
        *lane = v & 0xffff_ffff;
        c = v >> 32;
    }
    t[24] += c;

    // diagonal terms
    let mut c = 0u64;
    for i in 0..12 {
        let d = a32[i] * a32[i];
        let v = t[2 * i] + (d & 0xffff_ffff) + c;
        t[2 * i] = v & 0xffff_ffff;
        let v2 = t[2 * i + 1] + (d >> 32) + (v >> 32);
        t[2 * i + 1] = v2 & 0xffff_ffff;
        c = v2 >> 32;
    }
    t[24] += c;

    // Montgomery reduction, SOS shape
    for i in 0..12 {
        let m = t[i].wrapping_mul(INV32) & 0xffff_ffff;
        let mut carry = 0u64;
        for j in 0..12 {
            let (lo, hi) = mac32(t[i + j], m, p32[j], carry);
            t[i + j] = lo;
            carry = hi;
        }
        let mut k = i + 12;
        let mut c = carry;
        while c != 0 && k < 24 {
            let v = t[k] + c;
            t[k] = v & 0xffff_ffff;
            c = v >> 32;
            k += 1;
        }
        t[24] += c;
    }

    let mut r = [0u64; 6];
    for i in 0..6 {
        r[i] = t[12 + i * 2] | (t[12 + i * 2 + 1] << 32);
    }
    if t[24] != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

pub(crate) fn to_mont(a: &Fp) -> Fp {
    mont_mul(a, &R2)
}

pub(crate) fn from_mont(a: &Fp) -> Fp {
    mont_mul(a, &ONE)
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

fn modexp_bytes(base: &[u8], exp: &[u8]) -> Result<[u8; 48], ProgramError> {
    let modulus = limbs_to_be(&MODULUS);
    let params = BigModExpParams {
        base: base.as_ptr(),
        base_len: base.len() as u64,
        exponent: exp.as_ptr(),
        exponent_len: exp.len() as u64,
        modulus: modulus.as_ptr(),
        modulus_len: 48,
    };
    let mut out = [0u8; 48];
    let rc = unsafe {
        sys::sol_big_mod_exp(
            &params as *const BigModExpParams as *const u8,
            out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(out)
}

/// Modular exponentiation of a canonical-form element, returning canonical.
fn modexp(base: &Fp, exp: &[u8; 48]) -> Result<Fp, ProgramError> {
    Ok(be_to_limbs(&modexp_bytes(&limbs_to_be(base), exp)?))
}

/// Inverse of a Montgomery-form element, returned in Montgomery form.
fn inverse_mont(a: &Fp, exp_inv: &[u8; 48]) -> Result<Fp, ProgramError> {
    Ok(to_mont(&modexp(&from_mont(a), exp_inv)?))
}

fn expand_message_xmd(msg: &[u8]) -> [[u8; 32]; 4] {
    use solana_program::hash::hashv;

    let z_pad = [0u8; 64];
    let l_i_b = [0u8, 128];
    let dst_len = [DST_G1.len() as u8];

    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], DST_G1, &dst_len]).to_bytes();

    let mut blocks = [[0u8; 32]; 4];
    blocks[0] = hashv(&[&b0, &[1u8], DST_G1, &dst_len]).to_bytes();
    for i in 1..4 {
        let mut x = [0u8; 32];
        for j in 0..32 {
            x[j] = b0[j] ^ blocks[i - 1][j];
        }
        blocks[i] = hashv(&[&x, &[i as u8 + 1], DST_G1, &dst_len]).to_bytes();
    }
    blocks
}

/// hash_to_field for two Fp elements, canonical form.
fn hash_to_field(msg: &[u8]) -> Result<[Fp; 2], ProgramError> {
    let blocks = expand_message_xmd(msg);
    let one = [1u8];

    let mut wide = [0u8; 64];
    wide[..32].copy_from_slice(&blocks[0]);
    wide[32..].copy_from_slice(&blocks[1]);
    let u0 = be_to_limbs(&modexp_bytes(&wide, &one)?);

    wide[..32].copy_from_slice(&blocks[2]);
    wide[32..].copy_from_slice(&blocks[3]);
    let u1 = be_to_limbs(&modexp_bytes(&wide, &one)?);

    Ok([u0, u1])
}

struct Exps {
    inv: [u8; 48],
    legendre: [u8; 48],
    sqrt: [u8; 48],
}

impl Exps {
    fn new() -> Self {
        Self {
            inv: exp_inverse(),
            legendre: exp_legendre(),
            sqrt: exp_sqrt(),
        }
    }
}

/// Affine point on the 11-isogenous curve E', Montgomery form coordinates.
struct PointPrime {
    x: Fp,
    y: Fp,
}

fn is_one_canonical(a: &Fp) -> bool {
    a[0] == 1 && a[1..].iter().all(|&l| l == 0)
}

/// xi = 11: multiply by the SSWU non-residue with an addition chain.
fn mul_by_xi(a: &Fp) -> Fp {
    let a2 = add_mod(a, a);
    let a4 = add_mod(&a2, &a2);
    let a8 = add_mod(&a4, &a4);
    add_mod(&add_mod(&a8, &a2), a)
}

/// Simplified SWU map onto E', per RFC 9380 section 6.6.2 (variable time).
fn map_to_curve_sswu(u: &Fp, c_neg_b_over_a: &Fp, exps: &Exps) -> Result<PointPrime, ProgramError> {
    let um = to_mont(u);

    let usq = mont_sqr(&um);
    let xi_usq = mul_by_xi(&usq);
    let tv2 = add_mod(&mont_sqr(&xi_usq), &xi_usq);

    // tv2 == 0 has probability ~2^-381; the exceptional-case branch is omitted
    if is_zero(&tv2) {
        return Err(ProgramError::InvalidInstructionData);
    }

    let inv_tv2 = inverse_mont(&tv2, &exps.inv)?;
    let x1 = mont_mul(c_neg_b_over_a, &add_mod(&R, &inv_tv2));

    let gx = |x: &Fp| -> Fp {
        let xsq = mont_sqr(x);
        let x3 = mont_mul(&xsq, x);
        add_mod(&add_mod(&x3, &mont_mul(&SSWU_ELLP_A, x)), &SSWU_ELLP_B)
    };

    let gx1 = gx(&x1);
    let legendre = modexp(&from_mont(&gx1), &exps.legendre)?;

    let (x, gx_val) = if is_one_canonical(&legendre) {
        (x1, gx1)
    } else {
        let x2 = mont_mul(&xi_usq, &x1);
        let gx2 = gx(&x2);
        (x2, gx2)
    };

    let mut y = modexp(&from_mont(&gx_val), &exps.sqrt)?;

    // sgn0 correction: parity of y must match parity of u
    if (y[0] & 1) != (u[0] & 1) {
        y = neg_mod(&y);
    }

    Ok(PointPrime { x, y: to_mont(&y) })
}

/// Affine addition on E' (variable time; errors on the infinity outcome).
fn add_prime(p: &PointPrime, q: &PointPrime, exps: &Exps) -> Result<PointPrime, ProgramError> {
    let lambda = if p.x == q.x {
        if p.y != q.y || is_zero(&p.y) {
            return Err(ProgramError::InvalidInstructionData);
        }
        // doubling: (3x^2 + A) / 2y
        let xsq = mont_sqr(&p.x);
        let num = add_mod(&add_mod(&add_mod(&xsq, &xsq), &xsq), &SSWU_ELLP_A);
        let den = add_mod(&p.y, &p.y);
        mont_mul(&num, &inverse_mont(&den, &exps.inv)?)
    } else {
        let num = sub_mod(&q.y, &p.y);
        let den = sub_mod(&q.x, &p.x);
        mont_mul(&num, &inverse_mont(&den, &exps.inv)?)
    };

    let x3 = sub_mod(&sub_mod(&mont_sqr(&lambda), &p.x), &q.x);
    let y3 = sub_mod(&mont_mul(&lambda, &sub_mod(&p.x, &x3)), &p.y);
    Ok(PointPrime { x: x3, y: y3 })
}

/// Evaluate the four iso-11 polynomials with Knuth-adapted constants:
/// 27 multiplications against 51 for shared-nothing Horner. 
pub(crate) fn iso11_adapted(x: &Fp) -> (Fp, Fp, Fp, Fp) {
    let w = mont_mul(x, x);
    let mut xnum = add_mod(
        &mont_mul(&add_mod(x, &ISO11A_XNUM[0]), &add_mod(&w, &ISO11A_XNUM[1])),
        &ISO11A_XNUM[2],
    );
    let mut t = add_mod(&w, &ISO11A_XNUM[3]);
    xnum = add_mod(&mont_mul(&xnum, &t), &ISO11A_XNUM[4]);
    let mut t = add_mod(&w, &ISO11A_XNUM[5]);
    t = add_mod(&t, x);
    xnum = add_mod(&mont_mul(&xnum, &t), &ISO11A_XNUM[6]);
    let mut t = add_mod(&w, &ISO11A_XNUM[7]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    xnum = add_mod(&mont_mul(&xnum, &t), &ISO11A_XNUM[8]);
    let mut t = add_mod(&w, &ISO11A_XNUM[9]);
    xnum = add_mod(&mont_mul(&xnum, &t), &ISO11A_XNUM[10]);
    let xnum = mont_mul(&xnum, &ISO11A_XNUM[11]);

    let t = add_mod(&add_mod(&w, &mont_mul(&ISO11A_XDEN[0], x)), &ISO11A_XDEN[1]);
    let mut xden = add_mod(&mont_mul(&t, &add_mod(&w, &ISO11A_XDEN[2])), &ISO11A_XDEN[3]);
    let mut t = add_mod(&w, &ISO11A_XDEN[4]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    xden = add_mod(&mont_mul(&xden, &t), &ISO11A_XDEN[5]);
    let mut t = add_mod(&w, &ISO11A_XDEN[6]);
    t = add_mod(&t, x);
    xden = add_mod(&mont_mul(&xden, &t), &ISO11A_XDEN[7]);
    let mut t = add_mod(&w, &ISO11A_XDEN[8]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    xden = add_mod(&mont_mul(&xden, &t), &ISO11A_XDEN[9]);

    let mut ynum = add_mod(
        &mont_mul(&add_mod(x, &ISO11A_YNUM[0]), &add_mod(&w, &ISO11A_YNUM[1])),
        &ISO11A_YNUM[2],
    );
    let mut t = add_mod(&w, &ISO11A_YNUM[3]);
    t = add_mod(&t, x);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[4]);
    let mut t = add_mod(&w, &ISO11A_YNUM[5]);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[6]);
    let mut t = add_mod(&w, &ISO11A_YNUM[7]);
    t = add_mod(&t, x);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[8]);
    let mut t = add_mod(&w, &ISO11A_YNUM[9]);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[10]);
    let mut t = add_mod(&w, &ISO11A_YNUM[11]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[12]);
    let mut t = add_mod(&w, &ISO11A_YNUM[13]);
    ynum = add_mod(&mont_mul(&ynum, &t), &ISO11A_YNUM[14]);
    let ynum = mont_mul(&ynum, &ISO11A_YNUM[15]);

    let mut yden = add_mod(
        &mont_mul(&add_mod(x, &ISO11A_YDEN[0]), &add_mod(&w, &ISO11A_YDEN[1])),
        &ISO11A_YDEN[2],
    );
    let mut t = add_mod(&w, &ISO11A_YDEN[3]);
    t = add_mod(&t, x);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[4]);
    let mut t = add_mod(&w, &ISO11A_YDEN[5]);
    t = add_mod(&t, x);
    t = add_mod(&t, x);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[6]);
    let mut t = add_mod(&w, &ISO11A_YDEN[7]);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[8]);
    let mut t = add_mod(&w, &ISO11A_YDEN[9]);
    t = add_mod(&t, x);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[10]);
    let mut t = add_mod(&w, &ISO11A_YDEN[11]);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[12]);
    let mut t = add_mod(&w, &ISO11A_YDEN[13]);
    yden = add_mod(&mont_mul(&yden, &t), &ISO11A_YDEN[14]);

    (xnum, xden, ynum, yden)
}

fn horner(coeffs: &[[u64; 6]], x: &Fp) -> Fp {
    let mut acc = coeffs[coeffs.len() - 1];
    for c in coeffs[..coeffs.len() - 1].iter().rev() {
        acc = add_mod(&mont_mul(&acc, x), c);
    }
    acc
}

/// 11-isogeny from E' to E, affine in and out (Montgomery form).
fn iso_map(p: &PointPrime, exps: &Exps) -> Result<([u8; 48], [u8; 48]), ProgramError> {
    let (x_num, x_den, y_num, y_den) = iso11_adapted(&p.x);

    // batch inversion of both denominators
    let t = mont_mul(&x_den, &y_den);
    let t_inv = inverse_mont(&t, &exps.inv)?;
    let x_den_inv = mont_mul(&t_inv, &y_den);
    let y_den_inv = mont_mul(&t_inv, &x_den);

    let x = mont_mul(&x_num, &x_den_inv);
    let y = mont_mul(&p.y, &mont_mul(&y_num, &y_den_inv));

    Ok((limbs_to_be(&from_mont(&x)), limbs_to_be(&from_mont(&y))))
}

const BLS12_381_G1_BE: u64 = 5 | 0x80;
const OP_ADD: u64 = 0;

fn g1_add(a: &[u8; 96], b: &[u8; 96]) -> Result<[u8; 96], ProgramError> {
    let mut out = [0u8; 96];
    let rc = unsafe {
        sys::sol_curve_group_op(
            BLS12_381_G1_BE,
            OP_ADD,
            a.as_ptr(),
            b.as_ptr(),
            out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(out)
}

/// Multiplies by the effective cofactor with double-and-add over the g1 add
/// syscall, which skips the subgroup check that blocks the mul syscall here.
pub(crate) fn clear_cofactor(p: &[u8; 96]) -> Result<[u8; 96], ProgramError> {
    let mut acc = *p;
    for bit in (0..63).rev() {
        acc = g1_add(&acc, &acc)?;
        if (H_EFF >> bit) & 1 == 1 {
            acc = g1_add(&acc, p)?;
        }
    }
    Ok(acc)
}

pub(crate) fn validate(p: &[u8; 96]) -> Result<(), ProgramError> {
    let mut out = 0u8;
    let rc = unsafe { sys::sol_curve_validate_point(BLS12_381_G1_BE, p.as_ptr(), &mut out) };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(())
}

pub(crate) fn point_bytes(x: &[u8; 48], y: &[u8; 48]) -> [u8; 96] {
    let mut out = [0u8; 96];
    out[..48].copy_from_slice(x);
    out[48..].copy_from_slice(y);
    out
}

fn c_neg_b_over_a(exps: &Exps) -> Result<Fp, ProgramError> {
    let a_inv = inverse_mont(&SSWU_ELLP_A, &exps.inv)?;
    Ok(neg_mod(&mont_mul(&SSWU_ELLP_B, &a_inv)))
}

/// Stages, cumulative: 0 = hash_to_field, 1 = + both SSWU maps,
/// 2 = + E' add + isogeny, 3 = full with cofactor clearing and validation.
pub fn run(stage: u8, msg: &[u8]) -> Result<Vec<u8>, ProgramError> {
    let u = hash_to_field(msg)?;
    if stage == 0 {
        let mut out = Vec::with_capacity(96);
        out.extend_from_slice(&limbs_to_be(&u[0]));
        out.extend_from_slice(&limbs_to_be(&u[1]));
        return Ok(out);
    }

    let exps = Exps::new();
    let c = SSWU_C1_NEG_B_OVER_A;
    let p0 = map_to_curve_sswu(&u[0], &c, &exps)?;
    let p1 = map_to_curve_sswu(&u[1], &c, &exps)?;
    if stage == 1 {
        return Ok(limbs_to_be(&from_mont(&p0.x)).to_vec());
    }

    let sum = add_prime(&p0, &p1, &exps)?;
    let (x, y) = iso_map(&sum, &exps)?;
    let uncleared = point_bytes(&x, &y);
    if stage == 2 {
        return Ok(uncleared.to_vec());
    }

    let cleared = clear_cofactor(&uncleared)?;
    validate(&cleared)?;
    Ok(cleared.to_vec())
}


/// Fully unrolled 32-bit CIOS: constant indices elide bounds checks,
/// modulus limbs are immediates, loop control disappears.
#[rustfmt::skip]
fn mont_mul_cios32_unrolled(a: &Fp, b: &Fp) -> Fp {
    const M: u64 = 0xffff_ffff;
    let a0 = a[0] & M; let a1 = a[0] >> 32;
    let b0 = b[0] & M; let b1 = b[0] >> 32;
    let a2 = a[1] & M; let a3 = a[1] >> 32;
    let b2 = b[1] & M; let b3 = b[1] >> 32;
    let a4 = a[2] & M; let a5 = a[2] >> 32;
    let b4 = b[2] & M; let b5 = b[2] >> 32;
    let a6 = a[3] & M; let a7 = a[3] >> 32;
    let b6 = b[3] & M; let b7 = b[3] >> 32;
    let a8 = a[4] & M; let a9 = a[4] >> 32;
    let b8 = b[4] & M; let b9 = b[4] >> 32;
    let a10 = a[5] & M; let a11 = a[5] >> 32;
    let b10 = b[5] & M; let b11 = b[5] >> 32;
    let mut t = [0u64; 14];
    // round 0
    let mut c = 0u64;
    let v = t[0] + a0 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a0 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a0 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a0 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a0 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a0 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a0 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a0 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a0 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a0 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a0 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a0 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 1
    let mut c = 0u64;
    let v = t[0] + a1 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a1 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a1 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a1 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a1 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a1 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a1 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a1 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a1 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a1 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a1 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a1 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 2
    let mut c = 0u64;
    let v = t[0] + a2 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a2 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a2 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a2 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a2 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a2 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a2 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a2 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a2 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a2 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a2 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a2 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 3
    let mut c = 0u64;
    let v = t[0] + a3 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a3 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a3 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a3 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a3 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a3 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a3 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a3 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a3 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a3 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a3 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a3 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 4
    let mut c = 0u64;
    let v = t[0] + a4 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a4 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a4 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a4 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a4 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a4 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a4 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a4 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a4 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a4 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a4 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a4 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 5
    let mut c = 0u64;
    let v = t[0] + a5 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a5 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a5 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a5 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a5 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a5 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a5 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a5 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a5 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a5 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a5 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a5 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 6
    let mut c = 0u64;
    let v = t[0] + a6 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a6 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a6 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a6 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a6 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a6 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a6 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a6 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a6 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a6 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a6 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a6 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 7
    let mut c = 0u64;
    let v = t[0] + a7 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a7 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a7 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a7 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a7 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a7 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a7 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a7 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a7 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a7 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a7 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a7 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 8
    let mut c = 0u64;
    let v = t[0] + a8 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a8 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a8 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a8 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a8 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a8 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a8 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a8 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a8 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a8 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a8 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a8 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 9
    let mut c = 0u64;
    let v = t[0] + a9 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a9 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a9 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a9 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a9 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a9 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a9 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a9 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a9 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a9 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a9 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a9 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 10
    let mut c = 0u64;
    let v = t[0] + a10 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a10 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a10 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a10 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a10 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a10 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a10 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a10 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a10 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a10 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a10 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a10 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    // round 11
    let mut c = 0u64;
    let v = t[0] + a11 * b0 + c; t[0] = v & M; c = v >> 32;
    let v = t[1] + a11 * b1 + c; t[1] = v & M; c = v >> 32;
    let v = t[2] + a11 * b2 + c; t[2] = v & M; c = v >> 32;
    let v = t[3] + a11 * b3 + c; t[3] = v & M; c = v >> 32;
    let v = t[4] + a11 * b4 + c; t[4] = v & M; c = v >> 32;
    let v = t[5] + a11 * b5 + c; t[5] = v & M; c = v >> 32;
    let v = t[6] + a11 * b6 + c; t[6] = v & M; c = v >> 32;
    let v = t[7] + a11 * b7 + c; t[7] = v & M; c = v >> 32;
    let v = t[8] + a11 * b8 + c; t[8] = v & M; c = v >> 32;
    let v = t[9] + a11 * b9 + c; t[9] = v & M; c = v >> 32;
    let v = t[10] + a11 * b10 + c; t[10] = v & M; c = v >> 32;
    let v = t[11] + a11 * b11 + c; t[11] = v & M; c = v >> 32;
    let s = t[12] + c; t[12] = s & M; t[13] = s >> 32;
    let m = (t[0].wrapping_mul(0xfffcfffd)) & M;
    let v = t[0] + m * 0xffffaaab; let mut c = v >> 32;
    let v = t[1] + m * 0xb9feffff + c; t[0] = v & M; c = v >> 32;
    let v = t[2] + m * 0xb153ffff + c; t[1] = v & M; c = v >> 32;
    let v = t[3] + m * 0x1eabfffe + c; t[2] = v & M; c = v >> 32;
    let v = t[4] + m * 0xf6b0f624 + c; t[3] = v & M; c = v >> 32;
    let v = t[5] + m * 0x6730d2a0 + c; t[4] = v & M; c = v >> 32;
    let v = t[6] + m * 0xf38512bf + c; t[5] = v & M; c = v >> 32;
    let v = t[7] + m * 0x64774b84 + c; t[6] = v & M; c = v >> 32;
    let v = t[8] + m * 0x434bacd7 + c; t[7] = v & M; c = v >> 32;
    let v = t[9] + m * 0x4b1ba7b6 + c; t[8] = v & M; c = v >> 32;
    let v = t[10] + m * 0x397fe69a + c; t[9] = v & M; c = v >> 32;
    let v = t[11] + m * 0x1a0111ea + c; t[10] = v & M; c = v >> 32;
    let s = t[12] + c; t[11] = s & M; t[12] = t[13] + (s >> 32); t[13] = 0;
    let mut r = [0u64; 6];
    r[0] = t[0] | (t[1] << 32);
    r[1] = t[2] | (t[3] << 32);
    r[2] = t[4] | (t[5] << 32);
    r[3] = t[6] | (t[7] << 32);
    r[4] = t[8] | (t[9] << 32);
    r[5] = t[10] | (t[11] << 32);
    if t[12] != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

/// Product-scanning Montgomery multiply with 30-bit limbs and lazy
/// carries.
///
/// NOTE: reduces by 2^390, not 2^384.
const MASK30: u64 = (1 << 30) - 1;

fn to_limbs30(a: &Fp) -> [u64; 13] {
    let mut out = [0u64; 13];
    let mut bit = 0usize;
    for slot in out.iter_mut() {
        let word = bit / 64;
        let off = bit % 64;
        let mut v = a[word] >> off;
        if off > 34 && word + 1 < 6 {
            v |= a[word + 1] << (64 - off);
        }
        *slot = v & MASK30;
        bit += 30;
    }
    out
}

fn from_limbs30(a: &[u64; 13]) -> Fp {
    let mut out = [0u64; 6];
    for (i, &limb) in a.iter().enumerate() {
        let bit = i * 30;
        let word = bit / 64;
        let off = bit % 64;
        out[word] |= limb << off;
        if off > 34 && word + 1 < 6 {
            out[word + 1] |= limb >> (64 - off);
        }
    }
    out
}

fn mont_mul_comba30(a: &Fp, b: &Fp) -> Fp {
    let a30 = to_limbs30(a);
    let b30 = to_limbs30(b);
    let p30 = to_limbs30(&MODULUS);

    // full product, column scan: <= 13 products of < 2^60 per column
    // plus a < 2^35 running carry fits u64
    let mut t = [0u64; 26];
    let mut acc = 0u64;
    for col in 0usize..25 {
        let lo = col.saturating_sub(12);
        let hi = if col < 13 { col } else { 12 };
        for i in lo..=hi {
            acc += a30[i] * b30[col - i];
        }
        t[col] = acc & MASK30;
        acc >>= 30;
    }
    t[25] = acc;

    // reduction scan: fold 13 rounds of m*p column-wise
    let mut m = [0u64; 13];
    let mut acc = 0u64;
    for col in 0..13 {
        acc += t[col];
        for i in 0..col {
            acc += m[i] * p30[col - i];
        }
        m[col] = (acc & MASK30).wrapping_mul(INV30) & MASK30;
        acc += m[col] * p30[0];
        debug_assert_eq!(acc & MASK30, 0);
        acc >>= 30;
    }
    for col in 13..26 {
        if col < 26 {
            acc += t[col];
        }
        for i in (col - 12)..13 {
            acc += m[i] * p30[col - i];
        }
        t[col - 13] = acc & MASK30;
        acc >>= 30;
    }
    let mut r = from_limbs30(&t[..13].try_into().unwrap());
    if acc != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

/// -p^-1 mod 2^30 for the 30-bit-limb reduction.
const INV30: u64 = INV & MASK30;

/// p^2 in 32-bit lanes, the offset that keeps lazy Fp2 differences
/// non-negative.
pub(crate) const P2_WIDE: [u64; 24] = [
    0x1c718e39,
    0x26aa0000,
    0x76382eab,
    0x7ced6b1d,
    0x62113cfd,
    0x162c3383,
    0x3e71b743,
    0x66bf91ed,
    0x7091a049,
    0x292e85a8,
    0x86185c7b,
    0x1d68619c,
    0x0978ef01,
    0xf5314933,
    0x16ddca6e,
    0x50a62cfd,
    0x349e8bd0,
    0x66e59e49,
    0x0e7046b4,
    0xe2dc90e5,
    0xa22f25e9,
    0x4bd278ea,
    0xb8c35fc7,
    0x02a437a4,
];

pub(crate) const TWO_P2_WIDE: [u64; 24] = [
    0x38e31c72,
    0x4d540000,
    0xec705d56,
    0xf9dad63a,
    0xc42279fa,
    0x2c586706,
    0x7ce36e86,
    0xcd7f23da,
    0xe1234092,
    0x525d0b50,
    0x0c30b8f6,
    0x3ad0c339,
    0x12f1de02,
    0xea629266,
    0x2dbb94dd,
    0xa14c59fa,
    0x693d17a0,
    0xcdcb3c92,
    0x1ce08d68,
    0xc5b921ca,
    0x445e4bd3,
    0x97a4f1d5,
    0x7186bf8e,
    0x05486f49,
];

/// Full 768-bit product in 32-bit lanes, no reduction.
pub(crate) fn mul_wide32(a: &Fp, b: &Fp) -> [u64; 24] {
    let mut a32 = [0u64; 12];
    let mut b32 = [0u64; 12];
    for i in 0..6 {
        a32[i * 2] = a[i] & 0xffff_ffff;
        a32[i * 2 + 1] = a[i] >> 32;
        b32[i * 2] = b[i] & 0xffff_ffff;
        b32[i * 2 + 1] = b[i] >> 32;
    }
    let mut t = [0u64; 24];
    for i in 0..12 {
        let mut carry = 0u64;
        for j in 0..12 {
            let (lo, hi) = mac32(t[i + j], a32[i], b32[j], carry);
            t[i + j] = lo;
            carry = hi;
        }
        t[i + 12] = carry;
    }
    t
}

/// Montgomery reduction of a 768-bit lane vector; input below p * 2^384.
pub(crate) fn reduce_wide32(t: &mut [u64; 24]) -> Fp {
    let mut p32 = [0u64; 12];
    for i in 0..6 {
        p32[i * 2] = MODULUS[i] & 0xffff_ffff;
        p32[i * 2 + 1] = MODULUS[i] >> 32;
    }
    let mut carry2 = 0u64;
    for i in 0..12 {
        let m = t[i].wrapping_mul(INV32) & 0xffff_ffff;
        let mut carry = 0u64;
        for j in 0..12 {
            let (lo, hi) = mac32(t[i + j], m, p32[j], carry);
            t[i + j] = lo;
            carry = hi;
        }
        let v = t[i + 12] + carry2 + carry;
        t[i + 12] = v & 0xffff_ffff;
        carry2 = v >> 32;
    }
    let mut r = [0u64; 6];
    for i in 0..6 {
        r[i] = t[12 + i * 2] | (t[12 + i * 2 + 1] << 32);
    }
    if carry2 != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

pub(crate) fn wide_add(a: &[u64; 24], b: &[u64; 24]) -> [u64; 24] {
    let mut r = [0u64; 24];
    let mut c = 0u64;
    for i in 0..24 {
        let v = a[i] + b[i] + c;
        r[i] = v & 0xffff_ffff;
        c = v >> 32;
    }
    r
}

/// a - b for a >= b, lanes normalized.
pub(crate) fn wide_sub(a: &[u64; 24], b: &[u64; 24]) -> [u64; 24] {
    let mut r = [0u64; 24];
    let mut borrow = 0u64;
    for i in 0..24 {
        let v = a[i].wrapping_sub(b[i] + borrow);
        r[i] = v & 0xffff_ffff;
        borrow = (v >> 63) & 1;
    }
    r
}

/// CIOS variant: single interleaved pass, less memory traffic than SOS.
fn mont_mul_cios64(a: &Fp, b: &Fp) -> Fp {
    let mut t = [0u64; 8];
    for i in 0..6 {
        let mut carry = 0u64;
        for j in 0..6 {
            let (lo, hi) = mac(t[j], a[i], b[j], carry);
            t[j] = lo;
            carry = hi;
        }
        let (t6, c7) = adc(t[6], carry, 0);
        t[6] = t6;
        t[7] = c7;

        let m = t[0].wrapping_mul(INV);
        let (_, mut carry) = mac(t[0], m, MODULUS[0], 0);
        for j in 1..6 {
            let (lo, hi) = mac(t[j], m, MODULUS[j], carry);
            t[j - 1] = lo;
            carry = hi;
        }
        let (t5, c) = adc(t[6], carry, 0);
        t[5] = t5;
        t[6] = t[7] + c;
        t[7] = 0;
    }
    let mut r = [t[0], t[1], t[2], t[3], t[4], t[5]];
    if t[6] != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

const INV32: u64 = INV & 0xffff_ffff;

#[inline(always)]
fn mac32(acc: u64, a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = acc + a * b + carry;
    (t & 0xffff_ffff, t >> 32)
}

/// CIOS with 32-bit limbs: the multiply-accumulate needs no wide arithmetic.
fn mont_mul_cios32(a: &Fp, b: &Fp) -> Fp {
    let mut a32 = [0u64; 12];
    let mut b32 = [0u64; 12];
    let mut p32 = [0u64; 12];
    for i in 0..6 {
        a32[i * 2] = a[i] & 0xffff_ffff;
        a32[i * 2 + 1] = a[i] >> 32;
        b32[i * 2] = b[i] & 0xffff_ffff;
        b32[i * 2 + 1] = b[i] >> 32;
        p32[i * 2] = MODULUS[i] & 0xffff_ffff;
        p32[i * 2 + 1] = MODULUS[i] >> 32;
    }

    let mut t = [0u64; 14];
    for i in 0..12 {
        let ai = a32[i];
        let mut carry = 0u64;
        for j in 0..12 {
            let (lo, hi) = mac32(t[j], ai, b32[j], carry);
            t[j] = lo;
            carry = hi;
        }
        let s = t[12] + carry;
        t[12] = s & 0xffff_ffff;
        t[13] = s >> 32;

        let m = (t[0].wrapping_mul(INV32)) & 0xffff_ffff;
        let (_, mut carry) = mac32(t[0], m, p32[0], 0);
        for j in 1..12 {
            let (lo, hi) = mac32(t[j], m, p32[j], carry);
            t[j - 1] = lo;
            carry = hi;
        }
        let s = t[12] + carry;
        t[11] = s & 0xffff_ffff;
        t[12] = t[13] + (s >> 32);
        t[13] = 0;
    }

    let mut r = [0u64; 6];
    for i in 0..6 {
        r[i] = t[i * 2] | (t[i * 2 + 1] << 32);
    }
    if t[12] != 0 || geq(&r, &MODULUS) {
        r = sub_nocheck(&r, &MODULUS);
    }
    r
}

pub fn mul_bench(variant: u8, count: u64) -> u64 {
    let mul: fn(&Fp, &Fp) -> Fp = match variant {
        0 => mont_mul,
        1 => mont_mul_cios64,
        2 => mont_mul_cios32,
        3 => |a: &Fp, _: &Fp| mont_sqr_sos32(a),
        4 => |a: &Fp, _: &Fp| mont_mul(a, a),
        5 => mont_mul_comba30,
        _ => mont_mul_cios32_unrolled,
    };
    let mut acc = R2;
    let x = to_mont(&[3, 0, 0, 0, 0, 0]);
    for _ in 0..count {
        acc = mul(&core::hint::black_box(acc), &x);
    }
    acc[0]
}

pub fn mont_mul_bench(count: u64) -> u64 {
    let mut acc = R2;
    let x = to_mont(&[3, 0, 0, 0, 0, 0]);
    for _ in 0..count {
        acc = mont_mul(&core::hint::black_box(acc), &x);
    }
    acc[0]
}

// Witness-assisted variant: every result that big_mod_exp produced (inverses,
// square roots, branch selection) arrives as instruction data and gets checked
// with one or two multiplications. A wrong witness aborts; witnesses cannot
// steer the output point, which stays a pure function of the message.

const W_MAP: usize = 1 + 48 + 48;
const W_TOTAL: usize = 2 * W_MAP + 3 * 48;

use crate::g1_consts::SSWU_C1_NEG_B_OVER_A;

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

struct FieldElem {
    canonical: Fp,
    mont: Fp,
}

/// hash_to_field without modexp: split the 64-byte value at bit 256 and fold
/// with a precomputed 2^256 mod p.
fn hash_to_field_folded(msg: &[u8]) -> [FieldElem; 2] {
    let blocks = expand_message_xmd(msg);
    let mut out = [
        FieldElem { canonical: [0; 6], mont: [0; 6] },
        FieldElem { canonical: [0; 6], mont: [0; 6] },
    ];
    for (i, elem) in out.iter_mut().enumerate() {
        let mut hi = [0u8; 48];
        let mut lo = [0u8; 48];
        hi[16..].copy_from_slice(&blocks[i * 2]);
        lo[16..].copy_from_slice(&blocks[i * 2 + 1]);
        // canonical * Montgomery-form constant gives a canonical product
        let t = mont_mul(&be_to_limbs(&hi), &C256_MONT);
        let canonical = add_mod(&t, &be_to_limbs(&lo));
        elem.canonical = canonical;
        elem.mont = to_mont(&canonical);
    }
    out
}

fn gx_at(x: &Fp) -> Fp {
    let xsq = mont_sqr(x);
    let x3 = mont_mul(&xsq, x);
    add_mod(&add_mod(&x3, &mont_mul(&SSWU_ELLP_A, x)), &SSWU_ELLP_B)
}

/// The witness arrives in Montgomery form, so the check is one multiply.
pub(crate) fn check_inverse(v_m: &Fp, witness_m: &Fp) -> Result<Fp, ProgramError> {
    if mont_mul(v_m, witness_m) != R {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(*witness_m)
}

fn map_to_curve_witnessed(u: &FieldElem, wit: &[u8]) -> Result<PointPrime, ProgramError> {
    let flag = wit[0];
    let w_inv = wit48(&wit[1..49])?;
    let y_w = wit48(&wit[49..97])?;

    let usq = mont_sqr(&u.mont);
    let xi_usq = mul_by_xi(&usq);
    let tv2 = add_mod(&mont_sqr(&xi_usq), &xi_usq);
    if is_zero(&tv2) {
        return Err(ProgramError::InvalidInstructionData);
    }

    let inv_m = check_inverse(&tv2, &w_inv)?;
    let x1 = mont_mul(&SSWU_C1_NEG_B_OVER_A, &add_mod(&R, &inv_m));

    // gx2 = (Z u^2)^3 gx1 with Z a non-residue, so gx1 and gx2 always have
    // opposite quadratic characters: one sqrt witness proves its own branch.
    // A flag above 1 is non-canonical. The lone ambiguity is gx1 == 0 (which
    // forces gx2 == 0): both branches would then accept y = 0 and let the flag
    // steer the output between (x1,0) and (x2,0). blst takes x1 (is_square(0)
    // is true), so reject the x2 branch when gx == 0.
    if flag > 1 {
        return Err(ProgramError::InvalidInstructionData);
    }
    let (x, gx) = if flag == 0 {
        (x1, gx_at(&x1))
    } else {
        let x2 = mont_mul(&xi_usq, &x1);
        (x2, gx_at(&x2))
    };
    if flag == 1 && is_zero(&gx) {
        return Err(ProgramError::InvalidInstructionData);
    }

    let yw_m = to_mont(&y_w);
    if mont_sqr(&yw_m) != gx {
        return Err(ProgramError::InvalidInstructionData);
    }
    let mut y_canonical = y_w;

    if (y_canonical[0] & 1) != (u.canonical[0] & 1) {
        y_canonical = neg_mod(&y_canonical);
    }

    Ok(PointPrime { x, y: to_mont(&y_canonical) })
}

fn add_prime_witnessed(
    p: &PointPrime,
    q: &PointPrime,
    w_dx: &Fp,
) -> Result<PointPrime, ProgramError> {
    if p.x == q.x {
        return Err(ProgramError::InvalidInstructionData);
    }
    let dx = sub_mod(&q.x, &p.x);
    let inv_m = check_inverse(&dx, w_dx)?;
    let lambda = mont_mul(&sub_mod(&q.y, &p.y), &inv_m);
    let x3 = sub_mod(&sub_mod(&mont_sqr(&lambda), &p.x), &q.x);
    let y3 = sub_mod(&mont_mul(&lambda, &sub_mod(&p.x, &x3)), &p.y);
    Ok(PointPrime { x: x3, y: y3 })
}

fn iso_map_witnessed(
    p: &PointPrime,
    w_xd: &Fp,
    w_yd: &Fp,
) -> Result<([u8; 48], [u8; 48]), ProgramError> {
    let (x_num, x_den, y_num, y_den) = iso11_adapted(&p.x);

    let xd_inv = check_inverse(&x_den, w_xd)?;
    let yd_inv = check_inverse(&y_den, w_yd)?;

    let x = mont_mul(&x_num, &xd_inv);
    let y = mont_mul(&p.y, &mont_mul(&y_num, &yd_inv));
    Ok((limbs_to_be(&from_mont(&x)), limbs_to_be(&from_mont(&y))))
}


const DST_G1_NU: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_NU_POP_";

/// Single-element hash_to_field for the NU (encode_to_curve) suite.
fn hash_to_field_nu(msg: &[u8]) -> FieldElem {
    use solana_program::hash::hashv;
    let z_pad = [0u8; 64];
    let l_i_b = [0u8, 64];
    let dst_len = [DST_G1_NU.len() as u8];
    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], DST_G1_NU, &dst_len]).to_bytes();
    let b1 = hashv(&[&b0, &[1u8], DST_G1_NU, &dst_len]).to_bytes();
    let mut x = [0u8; 32];
    for j in 0..32 {
        x[j] = b0[j] ^ b1[j];
    }
    let b2 = hashv(&[&x, &[2u8], DST_G1_NU, &dst_len]).to_bytes();

    let mut hi = [0u8; 48];
    let mut lo = [0u8; 48];
    hi[16..].copy_from_slice(&b1);
    lo[16..].copy_from_slice(&b2);
    let t = mont_mul(&be_to_limbs(&hi), &C256_MONT);
    let canonical = add_mod(&t, &be_to_limbs(&lo));
    FieldElem { canonical, mont: to_mont(&canonical) }
}

/// Witnessed encode_to_curve (RFC 9380 NU): one map, no addition.
/// Blob: flag, w_inv, y, w_xd, w_yd then the message.
pub fn run_witnessed_nu(payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    const NU_TOTAL: usize = W_MAP + 2 * 48;
    if payload.len() < NU_TOTAL {
        return Err(ProgramError::InvalidInstructionData);
    }
    let (wits, msg) = payload.split_at(NU_TOTAL);
    let w_xd = wit48(&wits[W_MAP..W_MAP + 48])?;
    let w_yd = wit48(&wits[W_MAP + 48..])?;

    let u = hash_to_field_nu(msg);
    let p = map_to_curve_witnessed(&u, &wits[..W_MAP])?;
    let (x, y) = iso_map_witnessed(&p, &w_xd, &w_yd)?;
    let cleared = clear_cofactor(&point_bytes(&x, &y))?;
    validate(&cleared)?;
    Ok(cleared.to_vec())
}

pub fn run_witnessed(payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    if payload.len() < W_TOTAL {
        return Err(ProgramError::InvalidInstructionData);
    }
    let (wits, msg) = payload.split_at(W_TOTAL);

    let u = hash_to_field_folded(msg);
    let p0 = map_to_curve_witnessed(&u[0], &wits[..W_MAP])?;
    let p1 = map_to_curve_witnessed(&u[1], &wits[W_MAP..2 * W_MAP])?;

    let base = 2 * W_MAP;
    let w_dx = wit48(&wits[base..base + 48])?;
    let w_xd = wit48(&wits[base + 48..base + 96])?;
    let w_yd = wit48(&wits[base + 96..base + 144])?;

    let sum = add_prime_witnessed(&p0, &p1, &w_dx)?;
    let (x, y) = iso_map_witnessed(&sum, &w_xd, &w_yd)?;

    let cleared = clear_cofactor(&point_bytes(&x, &y))?;
    validate(&cleared)?;
    Ok(cleared.to_vec())
}

/// Host-side witness generation, mirroring the on-chain pipeline with the
/// expensive results computed via square-and-multiply.
#[cfg(not(target_os = "solana"))]
pub mod witness {
    use super::*;

    /// comba30 reduces by 2^390 where mont_mul reduces by 2^384, so
    /// comba30(a,b) * 2^6 must equal mont_mul(a,b) exactly.
    /// The adapted chains must agree with Horner over the original
    /// coefficient tables at arbitrary points.
    pub fn iso11_adapted_selftest() {
        let mut x = R2;
        for i in 0..50u64 {
            let expect = (
                horner(&ISO11_XNUM, &x),
                horner(&ISO11_XDEN, &x),
                horner(&ISO11_YNUM, &x),
                horner(&ISO11_YDEN, &x),
            );
            assert_eq!(iso11_adapted(&x), expect, "adapted iso11 diverged");
            x = add_mod(&mont_mul(&x, &R2), &[i, 1, 0, 0, 0, 0]);
        }
    }

    pub fn comba30_selftest() {
        let mut x: Fp = [1, 2, 3, 4, 5, 0];
        let mut y = R2;
        for _ in 0..100 {
            let mut lifted = mont_mul_comba30(&x, &y);
            for _ in 0..6 {
                lifted = add_mod(&lifted, &lifted);
            }
            assert_eq!(lifted, mont_mul(&x, &y), "comba30 disagrees with mont_mul");
            x = mont_mul(&x, &R2);
            y = add_mod(&mont_mul(&y, &y), &x);
        }
    }

    fn pow_mont(base: &Fp, exp_be: &[u8; 48]) -> Fp {
        let mut acc = R;
        for byte in exp_be {
            for bit in (0..8).rev() {
                acc = mont_mul(&acc, &acc);
                if (byte >> bit) & 1 == 1 {
                    acc = mont_mul(&acc, base);
                }
            }
        }
        acc
    }

    fn inverse(v_m: &Fp) -> [u8; 48] {
        limbs_to_be(&pow_mont(v_m, &exp_inverse()))
    }

    pub fn generate_nu(msg: &[u8]) -> Vec<u8> {
        let elem = hash_to_field_nu(msg);
        let (blob_map, point) = map_blob(&elem);
        let mut blob = blob_map;
        let x_den = horner(&ISO11_XDEN, &point.x);
        let y_den = horner(&ISO11_YDEN, &point.x);
        blob.extend_from_slice(&inverse(&x_den));
        blob.extend_from_slice(&inverse(&y_den));
        blob.extend_from_slice(msg);
        blob
    }

    /// One map's witness blob plus the mapped E' point.
    fn map_blob(elem: &FieldElem) -> (Vec<u8>, PointPrime) {
        let usq = mont_sqr(&elem.mont);
        let xi_usq = mul_by_xi(&usq);
        let tv2 = add_mod(&mont_sqr(&xi_usq), &xi_usq);
        assert!(!is_zero(&tv2));
        let w_inv = inverse(&tv2);
        let inv_m = be_to_limbs(&w_inv);
        let x1 = mont_mul(&SSWU_C1_NEG_B_OVER_A, &add_mod(&R, &inv_m));
        let gx1 = gx_at(&x1);
        let legendre = pow_mont(&gx1, &exp_legendre());
        let (flag, x, gx) = if is_zero(&gx1) || legendre == R {
            (0u8, x1, gx1)
        } else {
            let x2 = mont_mul(&xi_usq, &x1);
            (1u8, x2, gx_at(&x2))
        };
        let y_m = pow_mont(&gx, &exp_sqrt());
        assert_eq!(mont_sqr(&y_m), gx);
        let y_c = from_mont(&y_m);
        let mut y_final = y_c;
        if (y_final[0] & 1) != (elem.canonical[0] & 1) {
            y_final = neg_mod(&y_final);
        }
        let mut blob = vec![flag];
        blob.extend_from_slice(&w_inv);
        blob.extend_from_slice(&limbs_to_be(&y_c));
        (blob, PointPrime { x, y: to_mont(&y_final) })
    }

    pub fn generate(msg: &[u8]) -> Vec<u8> {
        let u = hash_to_field_folded(msg);
        let mut blob = Vec::with_capacity(W_TOTAL);
        let mut points = Vec::new();

        for elem in &u {
            let usq = mont_sqr(&elem.mont);
            let xi_usq = mul_by_xi(&usq);
            let tv2 = add_mod(&mont_sqr(&xi_usq), &xi_usq);
            assert!(!is_zero(&tv2));

            let w_inv = inverse(&tv2);
            let inv_m = be_to_limbs(&w_inv);
            let x1 = mont_mul(&SSWU_C1_NEG_B_OVER_A, &add_mod(&R, &inv_m));
            let gx1 = gx_at(&x1);

            // is_square(0) is true (matches blst / the on-chain zero guard).
            let legendre = pow_mont(&gx1, &exp_legendre());
            let (flag, x, gx) = if is_zero(&gx1) || legendre == R {
                (0u8, x1, gx1)
            } else {
                let x2 = mont_mul(&xi_usq, &x1);
                let gx2 = gx_at(&x2);
                (1u8, x2, gx2)
            };
            let y_m = pow_mont(&gx, &exp_sqrt());
            assert_eq!(mont_mul(&y_m, &y_m), gx);
            let y_c = from_mont(&y_m);

            let mut y_final = y_c;
            if (y_final[0] & 1) != (elem.canonical[0] & 1) {
                y_final = neg_mod(&y_final);
            }

            blob.push(flag);
            blob.extend_from_slice(&w_inv);
            blob.extend_from_slice(&limbs_to_be(&y_c));

            points.push(PointPrime { x, y: to_mont(&y_final) });
        }

        let dx = sub_mod(&points[1].x, &points[0].x);
        blob.extend_from_slice(&inverse(&dx));

        let inv_m = be_to_limbs(&inverse(&dx));
        let lambda = mont_mul(&sub_mod(&points[1].y, &points[0].y), &inv_m);
        let x3 = sub_mod(
            &sub_mod(&mont_mul(&lambda, &lambda), &points[0].x),
            &points[1].x,
        );
        let sum_x = x3;

        let x_den = horner(&ISO11_XDEN, &sum_x);
        let y_den = horner(&ISO11_YDEN, &sum_x);
        blob.extend_from_slice(&inverse(&x_den));
        blob.extend_from_slice(&inverse(&y_den));

        assert_eq!(blob.len(), W_TOTAL);
        blob
    }
}
