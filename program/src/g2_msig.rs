//! Witness-assisted RFC 9380 hash_to_G2 for BLS12-381 (min-pk suite
//! BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_).
//!
//! Same strategy as the G1 witness path but over Fp2: the transaction supplies
//! every inverse and square root as instruction data and the program verifies
//! each with a multiplication or two. Squareness of an Fp2 element is proven
//! through 1+i, a non-residue because 2 is a non-square mod p. Cofactor
//! clearing mirrors the Budroni-Pintore shape with the two mul-by-x chains
//! running through the g2 add syscall and the psi maps evaluated in-program.

use solana_program::program_error::ProgramError;

use crate::g1_consts::R;
use crate::g1_msig::{
    add_mod, be_to_limbs, from_mont, is_zero, limbs_to_be, mont_mul, neg_mod, sub_mod, sys,
    to_mont, wit48, Fp,
};
use crate::g2_consts::{
    C256_MONT, ISO3A_XDEN, ISO3A_XNUM, ISO3A_YDEN, ISO3A_YNUM, ISO3_XDEN, ISO3_YDEN,
    PSI2_X_C0, PSI_X_C1, PSI_Y,
    SSWU2_C1_NEG_B_OVER_A, SSWU2_ELLP_B,
};

const DST_G2: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
const BLS_X_ABS: u64 = 0xd201000000010000;

const BLS12_381_G2_BE: u64 = 6 | 0x80;
const OP_ADD: u64 = 0;
const OP_SUB: u64 = 1;
const POINT: usize = 192;

// blob: flag0, y0, flag1, y1, w_tv2_pair, w_dx, w_den_pair
const W_TOTAL: usize = 2 * (1 + 96) + 3 * 96;

#[derive(Clone, Copy, PartialEq, Debug)]
pub(crate) struct Fp2 {
    pub(crate) c0: Fp,
    pub(crate) c1: Fp,
}

pub(crate) const ZERO: Fp = [0; 6];
pub(crate) const ONE2: Fp2 = Fp2 { c0: R, c1: ZERO };

pub(crate) fn fp2(k: &[[u64; 6]; 2]) -> Fp2 {
    Fp2 { c0: k[0], c1: k[1] }
}

pub(crate) fn add2(a: &Fp2, b: &Fp2) -> Fp2 {
    Fp2 { c0: add_mod(&a.c0, &b.c0), c1: add_mod(&a.c1, &b.c1) }
}

pub(crate) fn sub2(a: &Fp2, b: &Fp2) -> Fp2 {
    Fp2 { c0: sub_mod(&a.c0, &b.c0), c1: sub_mod(&a.c1, &b.c1) }
}

pub(crate) fn neg2(a: &Fp2) -> Fp2 {
    Fp2 { c0: neg_mod(&a.c0), c1: neg_mod(&a.c1) }
}

pub(crate) fn is_zero2(a: &Fp2) -> bool {
    is_zero(&a.c0) && is_zero(&a.c1)
}

/// Karatsuba: works whenever the component products are valid mont_mul calls,
/// so also for canonical-times-Montgomery mixed-domain multiplication.
#[inline(always)]
pub(crate) fn mul2(a: &Fp2, b: &Fp2) -> Fp2 {
    let t0 = mont_mul(&a.c0, &b.c0);
    let t1 = mont_mul(&a.c1, &b.c1);
    let sa = add_mod(&a.c0, &a.c1);
    let sb = add_mod(&b.c0, &b.c1);
    let t2 = mont_mul(&sa, &sb);
    Fp2 {
        c0: sub_mod(&t0, &t1),
        c1: sub_mod(&sub_mod(&t2, &t0), &t1),
    }
}


#[inline(always)]
pub(crate) fn sq2(a: &Fp2) -> Fp2 {
    let s = add_mod(&a.c0, &a.c1);
    let d = sub_mod(&a.c0, &a.c1);
    let t = mont_mul(&a.c0, &a.c1);
    Fp2 {
        c0: mont_mul(&s, &d),
        c1: add_mod(&t, &t),
    }
}

pub(crate) fn to_mont2(a: &Fp2) -> Fp2 {
    Fp2 { c0: to_mont(&a.c0), c1: to_mont(&a.c1) }
}

pub(crate) fn from_mont2(a: &Fp2) -> Fp2 {
    Fp2 { c0: from_mont(&a.c0), c1: from_mont(&a.c1) }
}

pub(crate) fn wit96(bytes: &[u8]) -> Result<Fp2, ProgramError> {
    Ok(Fp2 {
        c0: wit48(&bytes[..48])?,
        c1: wit48(&bytes[48..96])?,
    })
}

pub(crate) fn sgn0_fp2(a: &Fp2) -> bool {
    // canonical form: sign of c0, falling back to c1 when c0 is zero
    let sign0 = a.c0[0] & 1 == 1;
    let zero0 = is_zero(&a.c0);
    let sign1 = a.c1[0] & 1 == 1;
    sign0 || (zero0 && sign1)
}

/// Affine point on the 3-isogenous curve E', Montgomery form.
pub(crate) struct Point2 {
    pub(crate) x: Fp2,
    pub(crate) y: Fp2,
}

fn expand_message_xmd_g2(msg: &[u8]) -> [[u8; 32]; 8] {
    use solana_program::hash::hashv;

    let z_pad = [0u8; 64];
    let l_i_b = [1u8, 0];
    let dst_len = [DST_G2.len() as u8];

    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], DST_G2, &dst_len]).to_bytes();

    let mut blocks = [[0u8; 32]; 8];
    blocks[0] = hashv(&[&b0, &[1u8], DST_G2, &dst_len]).to_bytes();
    for i in 1..8 {
        let mut x = [0u8; 32];
        for j in 0..32 {
            x[j] = b0[j] ^ blocks[i - 1][j];
        }
        blocks[i] = hashv(&[&x, &[i as u8 + 1], DST_G2, &dst_len]).to_bytes();
    }
    blocks
}

pub(crate) struct Elem2 {
    pub(crate) canonical: Fp2,
    pub(crate) mont: Fp2,
}

/// Reduce one 64-byte chunk: value = hi * 2^256 + lo with both halves < p.
pub(crate) fn fold(hi_block: &[u8; 32], lo_block: &[u8; 32]) -> (Fp, Fp) {
    let mut hi = [0u8; 48];
    let mut lo = [0u8; 48];
    hi[16..].copy_from_slice(hi_block);
    lo[16..].copy_from_slice(lo_block);
    // canonical * Montgomery-form constant gives a canonical product
    let t = mont_mul(&be_to_limbs(&hi), &C256_MONT);
    let canonical = add_mod(&t, &be_to_limbs(&lo));
    (canonical, to_mont(&canonical))
}

fn hash_to_field_g2(msg: &[u8]) -> [Elem2; 2] {
    let blocks = expand_message_xmd_g2(msg);
    let mut elems = [
        Elem2 { canonical: Fp2 { c0: ZERO, c1: ZERO }, mont: Fp2 { c0: ZERO, c1: ZERO } },
        Elem2 { canonical: Fp2 { c0: ZERO, c1: ZERO }, mont: Fp2 { c0: ZERO, c1: ZERO } },
    ];
    for (i, elem) in elems.iter_mut().enumerate() {
        let (c0, m0) = fold(&blocks[i * 4], &blocks[i * 4 + 1]);
        let (c1, m1) = fold(&blocks[i * 4 + 2], &blocks[i * 4 + 3]);
        elem.canonical = Fp2 { c0, c1 };
        elem.mont = Fp2 { c0: m0, c1: m1 };
    }
    elems
}

fn gx2_at(x: &Fp2) -> Fp2 {
    // x^3 + A x + B, with A of the form (0, a)
    let x3 = mul2(&sq2(x), x);
    let ax = mul_by_a2i(x);
    add2(&add2(&x3, &ax), &fp2(&SSWU2_ELLP_B))
}

/// The witness arrives in Montgomery form, so the check is one multiply.
pub(crate) fn check_inverse2(v: &Fp2, witness_m: &Fp2) -> Result<Fp2, ProgramError> {
    if mul2(v, witness_m) != ONE2 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(*witness_m)
}

/// xi2 = -(2 + i): (a + bi)(-2 - i) = (b - 2a) - (a + 2b) i. Adds only.
fn mul_by_xi2(v: &Fp2) -> Fp2 {
    let a2 = add_mod(&v.c0, &v.c0);
    let b2 = add_mod(&v.c1, &v.c1);
    Fp2 {
        c0: sub_mod(&v.c1, &a2),
        c1: neg_mod(&add_mod(&v.c0, &b2)),
    }
}

/// Multiply by 240 = 16 * (16 - 1) with an addition chain.
fn mul_fp_240(a: &Fp) -> Fp {
    let a2 = add_mod(a, a);
    let a4 = add_mod(&a2, &a2);
    let a8 = add_mod(&a4, &a4);
    let a15 = sub_mod(&add_mod(&a8, &a8), a);
    let t = add_mod(&a15, &a15);
    let t = add_mod(&t, &t);
    let t = add_mod(&t, &t);
    add_mod(&t, &t)
}

/// Multiply by A' = 240 i: (a + bi)(240 i) = -240 b + 240 a i.
fn mul_by_a2i(v: &Fp2) -> Fp2 {
    Fp2 {
        c0: neg_mod(&mul_fp_240(&v.c1)),
        c1: mul_fp_240(&v.c0),
    }
}

struct SswuPre {
    xi_usq: Fp2,
    tv2: Fp2,
}

fn sswu_pre(u: &Elem2) -> Result<SswuPre, ProgramError> {
    let usq = sq2(&u.mont);
    let xi_usq = mul_by_xi2(&usq);
    let tv2 = add2(&sq2(&xi_usq), &xi_usq);
    if is_zero2(&tv2) {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(SswuPre { xi_usq, tv2 })
}

/// Montgomery pair inversion: one witness w = (a*b)^-1 pins both inverses,
/// since a wrong w fails the single product check and inverses are unique.
fn check_pair_inverse(a: &Fp2, b: &Fp2, witness_m: &Fp2) -> Result<(Fp2, Fp2), ProgramError> {
    if mul2(&mul2(a, b), witness_m) != ONE2 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok((mul2(witness_m, b), mul2(witness_m, a)))
}

fn sswu_finish(
    u: &Elem2,
    pre: &SswuPre,
    inv: &Fp2,
    flag: u8,
    y_w: &Fp2,
) -> Result<Point2, ProgramError> {
    let x1 = mul2(&fp2(&SSWU2_C1_NEG_B_OVER_A), &add2(&ONE2, inv));

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
        (x1, gx2_at(&x1))
    } else {
        let x2 = mul2(&pre.xi_usq, &x1);
        (x2, gx2_at(&x2))
    };
    if flag == 1 && is_zero2(&gx) {
        return Err(ProgramError::InvalidInstructionData);
    }

    let yw_m = to_mont2(y_w);
    if sq2(&yw_m) != gx {
        return Err(ProgramError::InvalidInstructionData);
    }

    // Negation commutes with the Montgomery map, so flip the square root we
    // already converted rather than converting the canonical value a second time.
    let y = if sgn0_fp2(y_w) != sgn0_fp2(&u.canonical) {
        neg2(&yw_m)
    } else {
        yw_m
    };

    Ok(Point2 { x, y })
}

pub(crate) fn add_prime_witnessed(p: &Point2, q: &Point2, w_dx: &Fp2) -> Result<Point2, ProgramError> {
    if p.x == q.x {
        return Err(ProgramError::InvalidInstructionData);
    }
    let dx = sub2(&q.x, &p.x);
    let inv = check_inverse2(&dx, w_dx)?;
    let lambda = mul2(&sub2(&q.y, &p.y), &inv);
    let x3 = sub2(&sub2(&sq2(&lambda), &p.x), &q.x);
    let y3 = sub2(&mul2(&lambda, &sub2(&p.x, &x3)), &p.y);
    Ok(Point2 { x: x3, y: y3 })
}

fn horner2(coeffs: &[[[u64; 6]; 2]], x: &Fp2) -> Fp2 {
    let mut acc = fp2(&coeffs[coeffs.len() - 1]);
    for c in coeffs[..coeffs.len() - 1].iter().rev() {
        acc = add2(&mul2(&acc, x), &fp2(c));
    }
    acc
}

/// Adapted iso-3 evaluation: cubics as (y + g)(w + q1) + q0, the degree-2
/// denominator via its tiny coefficients (12 - 12i, -72i) as adds. Five
/// Fp2 multiplications plus one squaring against eleven for Horner.
fn iso3_adapted(x: &Fp2) -> (Fp2, Fp2, Fp2, Fp2) {
    fn m12(a: &Fp) -> Fp {
        let a2 = add_mod(a, a);
        let a4 = add_mod(&a2, &a2);
        add_mod(&add_mod(&a4, &a4), &a4)
    }
    let w = sq2(x);

    let cubic = |k: &[[[u64; 6]; 2]]| {
        add2(&mul2(&add2(x, &fp2(&k[0])), &add2(&w, &fp2(&k[1]))), &fp2(&k[2]))
    };
    let x_num = mul2(&cubic(&ISO3A_XNUM), &fp2(&ISO3A_XNUM[3]));
    let y_num = mul2(&cubic(&ISO3A_YNUM), &fp2(&ISO3A_YNUM[3]));
    let y_den = cubic(&ISO3A_YDEN);
    // x_den = w + (12 - 12i) x + k0
    let k1x = Fp2 {
        c0: m12(&add_mod(&x.c0, &x.c1)),
        c1: m12(&sub_mod(&x.c1, &x.c0)),
    };
    let x_den = add2(&add2(&w, &k1x), &fp2(&ISO3A_XDEN[1]));
    (x_num, x_den, y_num, y_den)
}

fn iso_map_witnessed(p: &Point2, w_den: &Fp2) -> Result<[u8; POINT], ProgramError> {
    let (x_num, x_den, y_num, y_den) = iso3_adapted(&p.x);

    let (xd_inv, yd_inv) = check_pair_inverse(&x_den, &y_den, w_den)?;

    let x = mul2(&x_num, &xd_inv);
    let y = mul2(&p.y, &mul2(&y_num, &yd_inv));
    Ok(point_bytes(&from_mont2(&x), &from_mont2(&y)))
}

/// Zcash uncompressed layout: x.c1 || x.c0 || y.c1 || y.c0, big-endian.
pub(crate) fn point_bytes(x: &Fp2, y: &Fp2) -> [u8; POINT] {
    let mut out = [0u8; POINT];
    out[..48].copy_from_slice(&limbs_to_be(&x.c1));
    out[48..96].copy_from_slice(&limbs_to_be(&x.c0));
    out[96..144].copy_from_slice(&limbs_to_be(&y.c1));
    out[144..].copy_from_slice(&limbs_to_be(&y.c0));
    out
}

fn parse_point(bytes: &[u8; POINT]) -> (Fp2, Fp2) {
    let x = Fp2 {
        c1: be_to_limbs(bytes[..48].try_into().unwrap()),
        c0: be_to_limbs(bytes[48..96].try_into().unwrap()),
    };
    let y = Fp2 {
        c1: be_to_limbs(bytes[96..144].try_into().unwrap()),
        c0: be_to_limbs(bytes[144..].try_into().unwrap()),
    };
    (x, y)
}

fn g2_group_op(op: u64, a: &[u8; POINT], b: &[u8; POINT]) -> Result<[u8; POINT], ProgramError> {
    let mut out = [0u8; POINT];
    let rc = unsafe {
        sys::sol_curve_group_op(BLS12_381_G2_BE, op, a.as_ptr(), b.as_ptr(), out.as_mut_ptr())
    };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(out)
}

fn g2_add(a: &[u8; POINT], b: &[u8; POINT]) -> Result<[u8; POINT], ProgramError> {
    g2_group_op(OP_ADD, a, b)
}

// The sub syscall skips the subgroup check like add, so it is safe on the
// pre-cleared cofactor intermediates and saves an in-program negation each call.
fn g2_sub(a: &[u8; POINT], b: &[u8; POINT]) -> Result<[u8; POINT], ProgramError> {
    g2_group_op(OP_SUB, a, b)
}

pub(crate) fn g2_validate(p: &[u8; POINT]) -> Result<(), ProgramError> {
    let mut out = 0u8;
    let rc = unsafe { sys::sol_curve_validate_point(BLS12_381_G2_BE, p.as_ptr(), &mut out) };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(())
}

fn neg_point(p: &[u8; POINT]) -> [u8; POINT] {
    let (x, y) = parse_point(p);
    point_bytes(&x, &neg2(&y))
}

/// [|x|]Q then negate, matching mul_by_x with BLS_X negative.
fn mul_by_x(q: &[u8; POINT]) -> Result<[u8; POINT], ProgramError> {
    let mut acc = *q;
    for bit in (0..63).rev() {
        acc = g2_add(&acc, &acc)?;
        if (BLS_X_ABS >> bit) & 1 == 1 {
            acc = g2_add(&acc, q)?;
        }
    }
    Ok(neg_point(&acc))
}

/// psi(x, y) = (conj(x) * (0, k_x), conj(y) * k_y). Constants are Montgomery,
/// coordinates canonical, so the mixed-domain products come out canonical.
fn psi(p: &[u8; POINT]) -> [u8; POINT] {
    let (x, y) = parse_point(p);
    // (a + bi) -> conj -> * (0, k): c0 = b*k, c1 = a*k
    let x_out = Fp2 {
        c0: mont_mul(&x.c1, &PSI_X_C1),
        c1: mont_mul(&x.c0, &PSI_X_C1),
    };
    let y_conj = Fp2 { c0: y.c0, c1: neg_mod(&y.c1) };
    let y_out = mul2(&y_conj, &fp2(&PSI_Y));
    point_bytes(&x_out, &y_out)
}

/// psi2(x, y) = (x * (k, 0), -y).
fn psi2(p: &[u8; POINT]) -> [u8; POINT] {
    let (x, y) = parse_point(p);
    let x_out = Fp2 {
        c0: mont_mul(&x.c0, &PSI2_X_C0),
        c1: mont_mul(&x.c1, &PSI2_X_C0),
    };
    point_bytes(&x_out, &neg2(&y))
}

/// Budroni-Pintore cofactor clearing, term for term as in the reference:
/// psi^2(2P) + [x^2 - x - 1]P + [x - 1]psi(P).
pub(crate) fn clear_cofactor(p: &[u8; POINT]) -> Result<[u8; POINT], ProgramError> {
    let t1 = mul_by_x(p)?;
    let t2 = psi(p);
    let p2 = psi2(&g2_add(p, p)?);
    let s2 = mul_by_x(&g2_add(&t1, &t2)?)?;

    let mut r = g2_add(&p2, &s2)?;
    r = g2_sub(&r, &t1)?;
    r = g2_sub(&r, &t2)?;
    r = g2_sub(&r, p)?;
    Ok(r)
}


const DST_G2_NU: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_NU_POP_";

/// Single-element hash_to_field for the NU (encode_to_curve) suite.
fn hash_to_field_nu(msg: &[u8]) -> Elem2 {
    use solana_program::hash::hashv;
    let z_pad = [0u8; 64];
    let l_i_b = [0u8, 128];
    let dst_len = [DST_G2_NU.len() as u8];
    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], DST_G2_NU, &dst_len]).to_bytes();
    let mut blocks = [[0u8; 32]; 4];
    blocks[0] = hashv(&[&b0, &[1u8], DST_G2_NU, &dst_len]).to_bytes();
    for i in 1..4 {
        let mut x = [0u8; 32];
        for j in 0..32 {
            x[j] = b0[j] ^ blocks[i - 1][j];
        }
        blocks[i] = hashv(&[&x, &[i as u8 + 1], DST_G2_NU, &dst_len]).to_bytes();
    }
    let (c0, m0) = fold(&blocks[0], &blocks[1]);
    let (c1, m1) = fold(&blocks[2], &blocks[3]);
    Elem2 { canonical: Fp2 { c0, c1 }, mont: Fp2 { c0: m0, c1: m1 } }
}

/// Witnessed encode_to_curve (RFC 9380 NU): one map, no addition.
/// Blob: flag, y, w_tv2, w_den then the message.
pub fn run_witnessed_nu(payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    const NU_TOTAL: usize = 1 + 96 + 96 + 96;
    if payload.len() < NU_TOTAL {
        return Err(ProgramError::InvalidInstructionData);
    }
    let (wits, msg) = payload.split_at(NU_TOTAL);
    let flag = wits[0];
    let y = wit96(&wits[1..97])?;
    let w_tv2 = wit96(&wits[97..193])?;
    let w_den = wit96(&wits[193..])?;

    let u = hash_to_field_nu(msg);
    let pre = sswu_pre(&u)?;
    let inv = check_inverse2(&pre.tv2, &w_tv2)?;
    let p = sswu_finish(&u, &pre, &inv, flag, &y)?;
    let uncleared = iso_map_witnessed(&p, &w_den)?;
    let cleared = clear_cofactor(&uncleared)?;
    g2_validate(&cleared)?;
    Ok(cleared.to_vec())
}

pub fn run_witnessed(payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    if payload.len() < W_TOTAL {
        return Err(ProgramError::InvalidInstructionData);
    }
    let (wits, msg) = payload.split_at(W_TOTAL);

    let flag0 = wits[0];
    let y0 = wit96(&wits[1..97])?;
    let flag1 = wits[97];
    let y1 = wit96(&wits[98..194])?;
    let w_tv2 = wit96(&wits[194..290])?;
    let w_dx = wit96(&wits[290..386])?;
    let w_den = wit96(&wits[386..482])?;

    let u = hash_to_field_g2(msg);
    let pre0 = sswu_pre(&u[0])?;
    let pre1 = sswu_pre(&u[1])?;
    let (inv0, inv1) = check_pair_inverse(&pre0.tv2, &pre1.tv2, &w_tv2)?;
    let p0 = sswu_finish(&u[0], &pre0, &inv0, flag0, &y0)?;
    let p1 = sswu_finish(&u[1], &pre1, &inv1, flag1, &y1)?;

    let sum = add_prime_witnessed(&p0, &p1, &w_dx)?;
    let uncleared = iso_map_witnessed(&sum, &w_den)?;

    let cleared = clear_cofactor(&uncleared)?;
    g2_validate(&cleared)?;
    Ok(cleared.to_vec())
}

/// Host-side witness generation mirroring the on-chain pipeline.
#[cfg(not(target_os = "solana"))]
pub mod witness {
    use super::*;
    use crate::g1_msig::{add_carryless, exp_inverse, exp_legendre, exp_sqrt, shr1};

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

    fn inv_fp(a: &Fp) -> Fp {
        pow_mont(a, &exp_inverse())
    }

    fn is_square_fp(a: &Fp) -> bool {
        is_zero(a) || pow_mont(a, &exp_legendre()) == R
    }

    fn sqrt_fp(a: &Fp) -> Fp {
        let s = pow_mont(a, &exp_sqrt());
        assert_eq!(mont_mul(&s, &s), *a, "not a square in Fp");
        s
    }

    /// Halve a residue: works in any representation.
    fn half(a: &Fp) -> Fp {
        if a[0] & 1 == 0 {
            shr1(a)
        } else {
            shr1(&add_carryless(a))
        }
    }

    fn norm(z: &Fp2) -> Fp {
        add_mod(&mont_mul(&z.c0, &z.c0), &mont_mul(&z.c1, &z.c1))
    }

    pub(crate) fn inv2(z: &Fp2) -> Fp2 {
        let n_inv = inv_fp(&norm(z));
        Fp2 {
            c0: mont_mul(&z.c0, &n_inv),
            c1: mont_mul(&neg_mod(&z.c1), &n_inv),
        }
    }

    pub(crate) fn is_square2(z: &Fp2) -> bool {
        is_square_fp(&norm(z))
    }

    /// Square root in Fp2 via the norm trick; input must be a square.
    pub(crate) fn sqrt2(z: &Fp2) -> Fp2 {
        if is_zero(&z.c1) {
            if is_square_fp(&z.c0) {
                return Fp2 { c0: sqrt_fp(&z.c0), c1: ZERO };
            }
            return Fp2 { c0: ZERO, c1: sqrt_fp(&neg_mod(&z.c0)) };
        }
        let delta = sqrt_fp(&norm(z));
        let mut t = half(&add_mod(&z.c0, &delta));
        if !is_square_fp(&t) {
            t = half(&sub_mod(&z.c0, &delta));
        }
        let x = sqrt_fp(&t);
        let y = mont_mul(&z.c1, &inv_fp(&add_mod(&x, &x)));
        let s = Fp2 { c0: x, c1: y };
        assert_eq!(sq2(&s), *z, "fp2 sqrt failed");
        s
    }

    pub(crate) fn push_fp2(blob: &mut Vec<u8>, z: &Fp2) {
        let c = from_mont2(z);
        blob.extend_from_slice(&limbs_to_be(&c.c0));
        blob.extend_from_slice(&limbs_to_be(&c.c1));
    }

    /// Serialize a Montgomery-form witness as-is (inverse witnesses).
    pub(crate) fn push_fp2_mont(blob: &mut Vec<u8>, z: &Fp2) {
        blob.extend_from_slice(&limbs_to_be(&z.c0));
        blob.extend_from_slice(&limbs_to_be(&z.c1));
    }

    // The other square root of gx: an equally valid witness that the sign
    // correction must resolve to the same output point.
    pub fn flip_first_sqrt(blob: &[u8]) -> Vec<u8> {
        let y = wit96(&blob[1..97]).unwrap();
        let mut out = blob[..1].to_vec();
        push_fp2(&mut out, &to_mont2(&neg2(&y)));
        out.extend_from_slice(&blob[97..]);
        out
    }

    pub fn generate_nu(msg: &[u8]) -> Vec<u8> {
        let u = hash_to_field_nu(msg);
        let pre = sswu_pre(&u).unwrap();
        let w_tv2 = inv2(&pre.tv2);
        let x1 = mul2(&fp2(&SSWU2_C1_NEG_B_OVER_A), &add2(&ONE2, &w_tv2));
        let gx1 = gx2_at(&x1);
        let (flag, x, gx) = if is_square2(&gx1) {
            (0u8, x1, gx1)
        } else {
            let x2 = mul2(&pre.xi_usq, &x1);
            let g = gx2_at(&x2);
            (1u8, x2, g)
        };
        let y = sqrt2(&gx);
        let mut y_canonical = from_mont2(&y);
        if sgn0_fp2(&y_canonical) != sgn0_fp2(&u.canonical) {
            y_canonical = neg2(&y_canonical);
        }
        let point = Point2 { x, y: to_mont2(&y_canonical) };
        let x_den = horner2(&ISO3_XDEN, &point.x);
        let y_den = horner2(&ISO3_YDEN, &point.x);
        let mut blob = vec![flag];
        push_fp2(&mut blob, &y);
        push_fp2_mont(&mut blob, &w_tv2);
        push_fp2_mont(&mut blob, &inv2(&mul2(&x_den, &y_den)));
        blob.extend_from_slice(msg);
        blob
    }

    pub fn generate(msg: &[u8]) -> Vec<u8> {
        let u = hash_to_field_g2(msg);
        let pre = [sswu_pre(&u[0]).unwrap(), sswu_pre(&u[1]).unwrap()];

        let w_tv2 = inv2(&mul2(&pre[0].tv2, &pre[1].tv2));
        let invs = [mul2(&w_tv2, &pre[1].tv2), mul2(&w_tv2, &pre[0].tv2)];

        let mut flags = [0u8; 2];
        let mut ys = [Fp2 { c0: ZERO, c1: ZERO }; 2];
        let mut points = Vec::new();
        for i in 0..2 {
            let x1 = mul2(&fp2(&SSWU2_C1_NEG_B_OVER_A), &add2(&ONE2, &invs[i]));
            let gx1 = gx2_at(&x1);
            let (flag, x, gx) = if is_square2(&gx1) {
                (0u8, x1, gx1)
            } else {
                let x2 = mul2(&pre[i].xi_usq, &x1);
                let gx2 = gx2_at(&x2);
                (1u8, x2, gx2)
            };
            let y = sqrt2(&gx);

            let mut y_canonical = from_mont2(&y);
            if sgn0_fp2(&y_canonical) != sgn0_fp2(&u[i].canonical) {
                y_canonical = neg2(&y_canonical);
            }

            flags[i] = flag;
            ys[i] = y;
            points.push(Point2 { x, y: to_mont2(&y_canonical) });
        }

        let dx = sub2(&points[1].x, &points[0].x);
        let dx_inv = inv2(&dx);

        let lambda = mul2(&sub2(&points[1].y, &points[0].y), &dx_inv);
        let x3 = sub2(&sub2(&sq2(&lambda), &points[0].x), &points[1].x);

        let x_den = horner2(&ISO3_XDEN, &x3);
        let y_den = horner2(&ISO3_YDEN, &x3);
        let w_den = inv2(&mul2(&x_den, &y_den));

        let mut blob = Vec::with_capacity(W_TOTAL);
        blob.push(flags[0]);
        push_fp2(&mut blob, &ys[0]);
        blob.push(flags[1]);
        push_fp2(&mut blob, &ys[1]);
        push_fp2_mont(&mut blob, &w_tv2);
        push_fp2_mont(&mut blob, &dx_inv);
        push_fp2_mont(&mut blob, &w_den);

        assert_eq!(blob.len(), W_TOTAL);
        blob
    }
}
