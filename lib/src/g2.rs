//! Witness-assisted RFC 9380 hash_to_G2 for BLS12-381 (min-pk suite
//! BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_).
//!
//! Same strategy as the G1 witness path but over Fp2: the transaction supplies
//! every inverse and square root as instruction data and the program verifies
//! each with a multiplication or two. Squareness of an Fp2 element is proven
//! through 1+i, a non-residue because 2 is a non-square mod p. Cofactor
//! clearing mirrors the Budroni-Pintore shape with the two mul-by-x chains
//! running through the g2 add syscall and the psi maps evaluated in-program.

use solana_program_error::ProgramError;
use alloc::vec::Vec;

use crate::fp2::*;
use crate::fp::{split_witness,
    add_mod, add_unreduced, be_to_limbs, from_mont, inv_xgcd, limbs_to_be,
    mont_mul, mont_sqr, neg_mod, sub_mod, sys, to_mont, wit48, Fp,
};
use crate::g1::check_inverse;
use crate::consts_g2::{
    C256_MONT, ISO3A_XDEN, ISO3A_XNUM, ISO3A_YDEN, ISO3A_YNUM,
    PSI2_X_C0, PSI_X_C1, PSI_Y,
    SSWU2_C1_NEG_B_OVER_A, SSWU2_ELLP_B,
};

const BLS_X_ABS: u64 = 0xd201000000010000;

const BLS12_381_G2_BE: u64 = 6 | 0x80;
const OP_ADD: u64 = 0;
const OP_SUB: u64 = 1;
const POINT: usize = 192;

// blob: flag0, y0, flag1, y1, tv2 witness(es), w_lambda, w_x, w_y
//
// Where an inverse feeds exactly one product (the E' slope, the iso-3
// output coordinates) the result is witnessed instead and pinned with a
// single multiply. 
#[cfg(feature = "wide-witness")]
const W_TV2: usize = 2 * 96;
#[cfg(not(feature = "wide-witness"))]
const W_TV2: usize = 96;
const W_TOTAL: usize = 2 * (1 + 96) + W_TV2 + 3 * 96;

// Compile time layout guards: the parse offsets below assume this exact
// blob shape in both feature configurations.
const _: () = assert!(W_TOTAL == 578 || W_TOTAL == 674);
const _: () = assert!(194 + W_TV2 + 3 * 96 == W_TOTAL);

pub(crate) fn expand_message_xmd_g2(dst: &[u8], msg: &[u8]) -> [[u8; 32]; 8] {
    use solana_sha256_hasher::hashv;

    let z_pad = [0u8; 64];
    let l_i_b = [1u8, 0];
    let dst_len = [dst.len() as u8];

    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], dst, &dst_len]).to_bytes();

    let mut blocks = [[0u8; 32]; 8];
    blocks[0] = hashv(&[&b0, &[1u8], dst, &dst_len]).to_bytes();
    for i in 1..8 {
        let mut x = [0u8; 32];
        for j in 0..32 {
            x[j] = b0[j] ^ blocks[i - 1][j];
        }
        blocks[i] = hashv(&[&x, &[i as u8 + 1], dst, &dst_len]).to_bytes();
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

fn hash_to_field_g2(dst: &[u8], msg: &[u8]) -> [Elem2; 2] {
    let blocks = expand_message_xmd_g2(dst, msg);
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

pub(crate) fn gx2_at(x: &Fp2) -> Fp2 {
    gx2_at_with_sq(x, &sq2(x))
}

/// x^3 + A x + B (A of the form (0, a)) given xsq = x^2, so callers that
/// already hold the square skip recomputing it.
#[inline(always)]
pub(crate) fn gx2_at_with_sq(x: &Fp2, xsq: &Fp2) -> Fp2 {
    let x3 = mul2(xsq, x);
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

/// One pair inversion witness pins both tv2 inverses through a single
/// product check, four mul2 for the smaller blob
#[cfg(not(feature = "wide-witness"))]
fn check_tv2_inverses(a: &Fp2, b: &Fp2, wit: &[u8]) -> Result<(Fp2, Fp2), ProgramError> {
    let w = wit96(wit)?;
    if mul2(&mul2(a, b), &w) != ONE2 {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok((mul2(&w, b), mul2(&w, a)))
}

/// Two direct inverse witnesses, two mul2 for 96 more blob bytes
#[cfg(feature = "wide-witness")]
fn check_tv2_inverses(a: &Fp2, b: &Fp2, wit: &[u8]) -> Result<(Fp2, Fp2), ProgramError> {
    let inv0 = check_inverse2(a, &wit96(&wit[..96])?)?;
    let inv1 = check_inverse2(b, &wit96(&wit[96..])?)?;
    Ok((inv0, inv1))
}

/// xi2 = -(2 + i): (a + bi)(-2 - i) = (b - 2a) - (a + 2b) i. Adds only.
pub(crate) fn mul_by_xi2(v: &Fp2) -> Fp2 {
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
pub(crate) fn mul_by_a2i(v: &Fp2) -> Fp2 {
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

    // The sqrt witness arrives in Montgomery form (a redc to read its sign
    // is cheaper than a to_mont to square-check it); negation commutes with
    // the Montgomery map, so the sign fix applies to the form we keep.
    let y_c = from_mont2(y_w);
    if sq2(y_w) != gx {
        return Err(ProgramError::InvalidInstructionData);
    }
    let y = if sgn0_fp2(&y_c) != sgn0_fp2(&u.canonical) {
        neg2(y_w)
    } else {
        *y_w
    };

    Ok(Point2 { x, y })
}

/// The slope itself is witnessed (not 1/dx): lambda * dx == dy pins it with
/// one multiply, and dx != 0 makes it unique.
pub(crate) fn add_prime_witnessed(p: &Point2, q: &Point2, w_lambda: &Fp2) -> Result<Point2, ProgramError> {
    if p.x == q.x {
        return Err(ProgramError::InvalidInstructionData);
    }
    let dx = sub2(&q.x, &p.x);
    let dy = sub2(&q.y, &p.y);
    if mul2(w_lambda, &dx) != dy {
        return Err(ProgramError::InvalidInstructionData);
    }
    let x3 = sub2(&sub2(&sq2(w_lambda), &p.x), &q.x);
    let y3 = sub2(&mul2(w_lambda, &sub2(&p.x, &x3)), &p.y);
    Ok(Point2 { x: x3, y: y3 })
}

fn horner2(coeffs: &[[[u64; 6]; 2]], x: &Fp2) -> Fp2 {
    let mut acc = fp2(&coeffs[coeffs.len() - 1]);
    for c in coeffs[..coeffs.len() - 1].iter().rev() {
        acc = add2(&mul2(&acc, x), &fp2(c));
    }
    acc
}

/// Multiply by 12 with an addition chain.
#[inline(always)]
fn m12(a: &Fp) -> Fp {
    let a2 = add_mod(a, a);
    let a4 = add_mod(&a2, &a2);
    add_mod(&add_mod(&a4, &a4), &a4)
}

/// One adapted cubic (x + g)(x^2 + q1) + q0 at x with xsq = x^2; the sums
/// ride unreduced, every one a mul2 operand or checked via mul2.
#[inline(always)]
fn iso3_cubic(x: &Fp2, xsq: &Fp2, k: &[[[u64; 6]; 2]]) -> Fp2 {
    add2_unreduced(
        &mul2(&add2_unreduced(x, &fp2(&k[0])), &add2_unreduced(xsq, &fp2(&k[1]))),
        &fp2(&k[2]),
    )
}

/// The iso-3 numerator pair at x given xsq = x^2: the adapted cubics times
/// their leading coefficients.
#[inline(always)]
fn iso3_numerators(x: &Fp2, xsq: &Fp2) -> (Fp2, Fp2) {
    let x_num = mul2(&iso3_cubic(x, xsq, &ISO3A_XNUM), &fp2(&ISO3A_XNUM[3]));
    let y_num = mul2(&iso3_cubic(x, xsq, &ISO3A_YNUM), &fp2(&ISO3A_YNUM[3]));
    (x_num, y_num)
}

/// Adapted iso-3 evaluation: cubics as (y + g)(w + q1) + q0, the degree-2
/// denominator via its tiny coefficients (12 - 12i, -72i) as adds. Five
/// Fp2 multiplications plus one squaring against eleven for Horner.
pub(crate) fn iso3_adapted(x: &Fp2) -> (Fp2, Fp2, Fp2, Fp2) {
    let w = sq2(x);
    let (x_num, y_num) = iso3_numerators(x, &w);
    let y_den = iso3_cubic(x, &w, &ISO3A_YDEN);
    // x_den = w + (12 - 12i) x + k0
    let k1x = Fp2 {
        c0: m12(&add_mod(&x.c0, &x.c1)),
        c1: m12(&sub_mod(&x.c1, &x.c0)),
    };
    let x_den = add2_unreduced(&add2_unreduced(&w, &k1x), &fp2(&ISO3A_XDEN[1]));
    (x_num, x_den, y_num, y_den)
}

/// The output coordinates are witnessed (Montgomery form) instead of the
/// denominator inverses: x * x_den == x_num pins x with one multiply per
/// coordinate. The RFC 9380 iso-3 numerators and denominators are coprime,
/// so both vanishing at once is impossible; a zero denominator (point maps
/// to infinity) is rejected as before.
fn iso_map_witnessed(p: &Point2, w_x: &Fp2, w_y: &Fp2) -> Result<[u8; POINT], ProgramError> {
    let (x_num, x_den, y_num, y_den) = iso3_adapted(&p.x);

    if is_zero2(&x_den) || is_zero2(&y_den) {
        return Err(ProgramError::InvalidInstructionData);
    }
    if mul2(w_x, &x_den) != x_num {
        return Err(ProgramError::InvalidInstructionData);
    }
    if mul2(w_y, &y_den) != mul2(&y_num, &p.y) {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(point_bytes(&from_mont2(w_x), &from_mont2(w_y)))
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
    // The syscall fills the whole point on success, so skip the zero-init;
    // the pointer escapes and LLVM cannot drop the memset on its own.
    let mut out = core::mem::MaybeUninit::<[u8; POINT]>::uninit();
    let rc = unsafe {
        sys::sol_curve_group_op(
            BLS12_381_G2_BE,
            op,
            a.as_ptr(),
            b.as_ptr(),
            out.as_mut_ptr() as *mut u8,
        )
    };
    if rc != 0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    // SAFETY: rc == 0 means the syscall wrote all POINT bytes
    Ok(unsafe { out.assume_init() })
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



/// Single-element hash_to_field for the NU (encode_to_curve) suite.
fn hash_to_field_nu(dst: &[u8], msg: &[u8]) -> Elem2 {
    use solana_sha256_hasher::hashv;
    let z_pad = [0u8; 64];
    let l_i_b = [0u8, 128];
    let dst_len = [dst.len() as u8];
    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], dst, &dst_len]).to_bytes();
    let mut blocks = [[0u8; 32]; 4];
    blocks[0] = hashv(&[&b0, &[1u8], dst, &dst_len]).to_bytes();
    for i in 1..4 {
        let mut x = [0u8; 32];
        for j in 0..32 {
            x[j] = b0[j] ^ blocks[i - 1][j];
        }
        blocks[i] = hashv(&[&x, &[i as u8 + 1], dst, &dst_len]).to_bytes();
    }
    let (c0, m0) = fold(&blocks[0], &blocks[1]);
    let (c1, m1) = fold(&blocks[2], &blocks[3]);
    Elem2 { canonical: Fp2 { c0, c1 }, mont: Fp2 { c0: m0, c1: m1 } }
}

/// Witnessed encode_to_curve (RFC 9380 NU): one map, no addition.
/// Blob: flag, y, w_tv2, w_x, w_y then the message.
pub fn encode_to_g2(dst: &[u8], payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    const NU_TOTAL: usize = 1 + 96 + 96 + 2 * 96;
    let (wits, msg) = split_witness(payload, NU_TOTAL)?;
    let flag = wits[0];
    let y = wit96(&wits[1..97])?;
    let w_tv2 = wit96(&wits[97..193])?;
    let w_x = wit96(&wits[193..289])?;
    let w_y = wit96(&wits[289..])?;

    let u = hash_to_field_nu(dst, msg);
    let pre = sswu_pre(&u)?;
    let inv = check_inverse2(&pre.tv2, &w_tv2)?;
    let p = sswu_finish(&u, &pre, &inv, flag, &y)?;
    let uncleared = iso_map_witnessed(&p, &w_x, &w_y)?;
    let cleared = clear_cofactor(&uncleared)?;
    g2_validate(&cleared)?;
    Ok(cleared.to_vec())
}

pub fn hash_to_g2(dst: &[u8], payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    hash_to_g2_prefix(dst, 3, payload)
}

/// Cumulative stage prefixes of the pipeline for stage by stage CU
/// measurement; stage 3 is the full hash
#[doc(hidden)]
pub fn hash_to_g2_prefix(dst: &[u8], stage: u8, payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    let (wits, msg) = split_witness(payload, W_TOTAL)?;

    let flag0 = wits[0];
    let y0 = wit96(&wits[1..97])?;
    let flag1 = wits[97];
    let y1 = wit96(&wits[98..194])?;
    let tv2_witness = &wits[194..194 + W_TV2];
    let w_lambda = wit96(&wits[194 + W_TV2..290 + W_TV2])?;
    let w_x = wit96(&wits[290 + W_TV2..386 + W_TV2])?;
    let w_y = wit96(&wits[386 + W_TV2..482 + W_TV2])?;

    let u = hash_to_field_g2(dst, msg);
    if stage == 0 {
        return Ok(limbs_to_be(&u[0].canonical.c0).to_vec());
    }

    let pre0 = sswu_pre(&u[0])?;
    let pre1 = sswu_pre(&u[1])?;
    let (inv0, inv1) = check_tv2_inverses(&pre0.tv2, &pre1.tv2, tv2_witness)?;
    let p0 = sswu_finish(&u[0], &pre0, &inv0, flag0, &y0)?;
    let p1 = sswu_finish(&u[1], &pre1, &inv1, flag1, &y1)?;
    if stage == 1 {
        return Ok(limbs_to_be(&p1.x.c0).to_vec());
    }

    let sum = add_prime_witnessed(&p0, &p1, &w_lambda)?;
    let uncleared = iso_map_witnessed(&sum, &w_x, &w_y)?;
    if stage == 2 {
        return Ok(uncleared.to_vec());
    }

    let cleared = clear_cofactor(&uncleared)?;
    g2_validate(&cleared)?;
    Ok(cleared.to_vec())
}

// ---------------------------------------------------------------------------
// Compact witness: 145 bytes against the default blob's 578.
//
// blob: flags (two branch bits), c0(y0), c0(y1), w
//
// Two ideas shrink the blob. A square root over Fp2 travels as its real
// part alone: with gx = g0 + g1 i and y = c0 + c1 i, the imaginary half of
// y^2 == gx forces c1 = g1 / (2 c0), so shipping c0 and checking the real
// half (c0^2 - c1^2 == g0) pins the root. And every inverse the pipeline
// needs -- the two tv2 inverses, the two 2c0 divisors, and per map the
// iso-3 denominator pair -- collapses into one Montgomery batch: the blob
// carries a single Fp witness w = (e1 ... e8)^-1 whose product check pins
// every unpacked inverse. Fp2 inverses ride their norms (inv z =
// conj z / norm z, and norm z = 0 only at z = 0 since -1 is a non-residue),
// so the batch stays in Fp.
//
// The iso-3 denominators must join the batch before any inverse exists, so
// they are evaluated homogenized over the fractional SSWU candidate
// x' = n / tv2 (xden_h = tv2^2 xden(x'), yden_h = tv2^3 yden(x')); the
// numerators wait for the affine x'. Each map then lands on E on its own
// and the E-side addition rides the g2 add syscall instead of a slope
// witness, which also drops the equal-x abort of the witnessed E' add.
//
// Steering room: none. w is pinned by the product check; c0 by the
// real-part check, and both roots of gx collapse to one output through the
// sign rule as before; a lied branch flag either changes the batch elements
// out from under w or leaves the real-part check unsatisfiable (the
// wrong branch's gx is a non-square by the SSWU non-residue construction).
// A zero anywhere (u = 0, iso poles, a zero denominator) zeroes the batch
// product and aborts. Like the rest of the pipeline this trades the RFC's
// exceptional cases for aborts on measure-zero inputs; the one new case is
// gx square with a purely imaginary root (c0 = 0), which cannot be encoded.

const C_TOTAL: usize = 1 + 3 * 48;
// witness-free-inverse layout: the compact blob minus its trailing w
const C_XGCD_TOTAL: usize = 1 + 2 * 48;

pub(crate) fn norm_fp2(z: &Fp2) -> Fp {
    add_mod(&mont_sqr(&z.c0), &mont_sqr(&z.c1))
}

/// Replace every element with its inverse in place (all Montgomery form)
/// through one shared inversion. With a witness the product check pins it
/// (one multiply); without one the product is inverted in-program by
/// binary xgcd, trading the 48 witness bytes for shift-and-subtract CU.
/// Any zero element zeroes the product and fails either path.
/// inline(never): the prefix table lives in its own SBF frame.
#[inline(never)]
fn batch_inverse<const N: usize>(elems: &mut [Fp; N], w: Option<&Fp>) -> Result<(), ProgramError> {
    let mut prefix = [ZERO; N];
    prefix[0] = elems[0];
    for i in 1..N {
        prefix[i] = mont_mul(&prefix[i - 1], &elems[i]);
    }
    let w = match w {
        Some(w) => check_inverse(&prefix[N - 1], w)?,
        None => to_mont(&inv_xgcd(&from_mont(&prefix[N - 1]))?),
    };
    let mut acc = w;
    for i in (1..N).rev() {
        let e = elems[i];
        elems[i] = mont_mul(&acc, &prefix[i - 1]);
        acc = mont_mul(&acc, &e);
    }
    elems[0] = acc;
    Ok(())
}

/// Fp2 inverse assembled from the batched inverse of its norm. Negation
/// follows the multiply so unreduced components are legal inputs.
fn inv2_from_norm(z: &Fp2, norm_inv: &Fp) -> Fp2 {
    Fp2 {
        c0: mont_mul(&z.c0, norm_inv),
        c1: neg_mod(&mont_mul(&z.c1, norm_inv)),
    }
}

/// Per-map values both compact paths derive before any inverse exists: the
/// fractional SSWU candidate numerator (denominator tv2) and the iso-3
/// denominators evaluated homogenized over it (xden_h = tv2^2 xden(x'),
/// yden_h = tv2^3 yden(x')), so their norms can join the batch. Shared by
/// the verifier and the host generator: the batch layout is defined once.
struct MapParts {
    n: Fp2,
    dsq: Fp2,
    dcb: Fp2,
    xden_h: Fp2,
    yden_h: Fp2,
}

/// inline(never): the intermediates live in their own SBF frame.
#[inline(never)]
fn compact_map_parts(pre: &SswuPre, flag: u8) -> MapParts {
    let mut n = mul2(&fp2(&SSWU2_C1_NEG_B_OVER_A), &add2_unreduced(&pre.tv2, &ONE2));
    if flag == 1 {
        n = mul2(&pre.xi_usq, &n);
    }
    let d = &pre.tv2;
    let dsq = sq2(d);
    let dcb = mul2(&dsq, d);
    let nsq = sq2(&n);
    let nd = mul2(&n, d);
    // sums ride unreduced: they feed mont_sqr norms and inv2_from_norm
    // xden_h = n^2 + (12 - 12i) n d + k0 d^2
    let k1nd = Fp2 {
        c0: m12(&add_mod(&nd.c0, &nd.c1)),
        c1: m12(&sub_mod(&nd.c1, &nd.c0)),
    };
    let xden_h = add2_unreduced(&add2_unreduced(&nsq, &k1nd), &mul2(&fp2(&ISO3A_XDEN[1]), &dsq));
    // the adapted monic cubic (x + g)(x^2 + q1) + q0, homogenized
    let yden_h = add2_unreduced(
        &mul2(
            &add2_unreduced(&n, &mul2(&fp2(&ISO3A_YDEN[0]), d)),
            &add2_unreduced(&nsq, &mul2(&fp2(&ISO3A_YDEN[1]), &dsq)),
        ),
        &mul2(&fp2(&ISO3A_YDEN[2]), &dcb),
    );
    MapParts { n, dsq, dcb, xden_h, yden_h }
}

/// The per-map denominators after batch unpacking, bundled for the tail.
struct IsoDen<'a> {
    dsq: &'a Fp2,
    dcb: &'a Fp2,
    inv_xden: Fp2,
    inv_yden: Fp2,
}

/// One map's tail: recover the sqrt from its real part, fix the sign, and
/// land the point on E through the iso with pre-inverted denominators.
#[inline(never)]
fn sswu_iso_compact(
    u: &Elem2,
    x: &Fp2,
    flag: u8,
    c0: &Fp,
    inv_2c0: &Fp,
    den: &IsoDen,
) -> Result<[u8; POINT], ProgramError> {
    let xsq = sq2(x);
    let gx = gx2_at_with_sq(x, &xsq);
    // gx == 0 forces both branches to zero; blst takes x1, so the x2 branch
    // is rejected (its c0 = 0 case aborts in the batch before this)
    if flag == 1 && is_zero2(&gx) {
        return Err(ProgramError::InvalidInstructionData);
    }

    // c1 = g1 / (2 c0) satisfies the imaginary half of y^2 == gx by
    // construction; the real half, as a difference of squares, is the
    // remaining unique-answer check
    let c1 = mont_mul(&gx.c1, inv_2c0);
    if mont_mul(&add_unreduced(c0, &c1), &sub_mod(c0, &c1)) != gx.c0 {
        return Err(ProgramError::InvalidInstructionData);
    }
    let y_m = Fp2 { c0: *c0, c1 };
    // c0 != 0 (its double is a batch element), so sgn0(y) is c0's parity
    let y = if (from_mont(c0)[0] & 1 == 1) != sgn0_fp2(&u.canonical) {
        neg2(&y_m)
    } else {
        y_m
    };

    // numerators only: the denominators arrived through the batch
    let (x_num, y_num) = iso3_numerators(x, &xsq);
    let x_out = mul2(&mul2(&x_num, den.dsq), &den.inv_xden);
    let y_out = mul2(&mul2(&mul2(&y, &y_num), den.dcb), &den.inv_yden);
    Ok(point_bytes(&from_mont2(&x_out), &from_mont2(&y_out)))
}

/// hash_to_g2 against the 145-byte compact blob; same suite, same output
/// bytes, roughly the CU of the default path plus the traded inversions.
pub fn hash_to_g2_compact(dst: &[u8], payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    let (wits, msg) = split_witness(payload, C_TOTAL)?;
    let w = wit48(&wits[C_XGCD_TOTAL..C_TOTAL])?;
    hash_to_g2_compact_inner(dst, msg, wits, Some(&w))
}

/// hash_to_g2 against the 97-byte blob: flags and the two root halves
/// only. The batched inverse is computed in-program by binary xgcd
/// instead of witnessed, trading its 48 bytes for shift-and-subtract CU.
pub fn hash_to_g2_compact_xgcd(dst: &[u8], payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    let (wits, msg) = split_witness(payload, C_XGCD_TOTAL)?;
    hash_to_g2_compact_inner(dst, msg, wits, None)
}

fn hash_to_g2_compact_inner(
    dst: &[u8],
    msg: &[u8],
    wits: &[u8],
    w: Option<&Fp>,
) -> Result<Vec<u8>, ProgramError> {
    let flags = wits[0];
    if flags > 3 {
        return Err(ProgramError::InvalidInstructionData);
    }
    let c0s = [wit48(&wits[1..49])?, wit48(&wits[49..97])?];

    let u = hash_to_field_g2(dst, msg);
    let pre = [sswu_pre(&u[0])?, sswu_pre(&u[1])?];
    let parts = [
        compact_map_parts(&pre[0], flags & 1),
        compact_map_parts(&pre[1], (flags >> 1) & 1),
    ];

    // after the in-place batch these hold the inverses, element for element
    let mut inv = [
        add_unreduced(&c0s[0], &c0s[0]),
        add_unreduced(&c0s[1], &c0s[1]),
        norm_fp2(&pre[0].tv2),
        norm_fp2(&pre[1].tv2),
        norm_fp2(&parts[0].xden_h),
        norm_fp2(&parts[0].yden_h),
        norm_fp2(&parts[1].xden_h),
        norm_fp2(&parts[1].yden_h),
    ];
    batch_inverse(&mut inv, w)?;

    let mut points = [[0u8; POINT]; 2];
    for i in 0..2 {
        let x = mul2(&parts[i].n, &inv2_from_norm(&pre[i].tv2, &inv[2 + i]));
        let den = IsoDen {
            dsq: &parts[i].dsq,
            dcb: &parts[i].dcb,
            inv_xden: inv2_from_norm(&parts[i].xden_h, &inv[4 + 2 * i]),
            inv_yden: inv2_from_norm(&parts[i].yden_h, &inv[5 + 2 * i]),
        };
        points[i] = sswu_iso_compact(&u[i], &x, (flags >> i) & 1, &c0s[i], &inv[i], &den)?;
    }

    let sum = g2_add(&points[0], &points[1])?;
    let cleared = clear_cofactor(&sum)?;
    g2_validate(&cleared)?;
    Ok(cleared.to_vec())
}

/// Host-side witness generation mirroring the on-chain pipeline.
#[cfg(not(target_os = "solana"))]
pub mod witness {
    use super::*;
    use crate::consts_g1::R;
    use crate::consts_g2::{ISO3_XDEN, ISO3_XNUM, ISO3_YDEN, ISO3_YNUM};
    use crate::fp::{exp_inverse, exp_legendre, exp_sqrt, half_mod, is_zero};
    use alloc::vec;

    fn pow_mont(base: &Fp, exp_be: &[u8; 48]) -> Fp {
        let mut table = [R; 16];
        table[1] = *base;
        for i in 2..16 {
            table[i] = mont_mul(&table[i - 1], base);
        }
        let mut acc = R;
        for byte in exp_be {
            for nib in [byte >> 4, byte & 0xf] {
                for _ in 0..4 {
                    acc = mont_mul(&acc, &acc);
                }
                if nib != 0 {
                    acc = mont_mul(&acc, &table[nib as usize]);
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

    pub(crate) fn inv2(z: &Fp2) -> Fp2 {
        let n_inv = inv_fp(&norm_fp2(z));
        Fp2 {
            c0: mont_mul(&z.c0, &n_inv),
            c1: mont_mul(&neg_mod(&z.c1), &n_inv),
        }
    }

    pub(crate) fn is_square2(z: &Fp2) -> bool {
        is_square_fp(&norm_fp2(z))
    }

    /// Square root in Fp2 via the norm trick; input must be a square.
    pub(crate) fn sqrt2(z: &Fp2) -> Fp2 {
        if is_zero(&z.c1) {
            if is_square_fp(&z.c0) {
                return Fp2 { c0: sqrt_fp(&z.c0), c1: ZERO };
            }
            return Fp2 { c0: ZERO, c1: sqrt_fp(&neg_mod(&z.c0)) };
        }
        let delta = sqrt_fp(&norm_fp2(z));
        let mut t = half_mod(&add_mod(&z.c0, &delta));
        if !is_square_fp(&t) {
            t = half_mod(&sub_mod(&z.c0, &delta));
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
        push_fp2_mont(&mut out, &neg2(&y));
        out.extend_from_slice(&blob[97..]);
        out
    }

    /// SSWU branch selection at the affine candidate: x1 when gx1 is a
    /// square, else x2 = xi u^2 x1.
    fn select_sswu_branch(pre: &SswuPre, inv_tv2: &Fp2) -> (u8, Fp2, Fp2) {
        let x1 = mul2(&fp2(&SSWU2_C1_NEG_B_OVER_A), &add2(&ONE2, inv_tv2));
        let gx1 = gx2_at(&x1);
        if is_square2(&gx1) {
            (0, x1, gx1)
        } else {
            let x2 = mul2(&pre.xi_usq, &x1);
            let gx2 = gx2_at(&x2);
            (1, x2, gx2)
        }
    }

    pub fn generate_nu(msg: &[u8]) -> Vec<u8> {
        let u = hash_to_field_nu(crate::dst::G2_NU, msg);
        let pre = sswu_pre(&u).unwrap();
        let w_tv2 = inv2(&pre.tv2);
        let (flag, x, gx) = select_sswu_branch(&pre, &w_tv2);
        let y = sqrt2(&gx);
        let mut y_canonical = from_mont2(&y);
        if sgn0_fp2(&y_canonical) != sgn0_fp2(&u.canonical) {
            y_canonical = neg2(&y_canonical);
        }
        let point = Point2 { x, y: to_mont2(&y_canonical) };
        let (x_out, y_out) = iso_outputs(&point);
        let mut blob = vec![flag];
        push_fp2_mont(&mut blob, &y);
        push_fp2_mont(&mut blob, &w_tv2);
        push_fp2_mont(&mut blob, &x_out);
        push_fp2_mont(&mut blob, &y_out);
        blob.extend_from_slice(msg);
        blob
    }

    /// The iso-3 image of an E' point, Montgomery form: the witnessed
    /// output coordinates.
    fn iso_outputs(p: &Point2) -> (Fp2, Fp2) {
        let x_num = horner2(&ISO3_XNUM, &p.x);
        let y_num = horner2(&ISO3_YNUM, &p.x);
        let x_den = horner2(&ISO3_XDEN, &p.x);
        let y_den = horner2(&ISO3_YDEN, &p.x);
        let w = inv2(&mul2(&x_den, &y_den));
        let x_out = mul2(&x_num, &mul2(&w, &y_den));
        let y_out = mul2(&mul2(&p.y, &y_num), &mul2(&w, &x_den));
        (x_out, y_out)
    }

    pub fn generate(msg: &[u8]) -> Vec<u8> {
        let u = hash_to_field_g2(crate::dst::G2_RO, msg);
        let pre = [sswu_pre(&u[0]).unwrap(), sswu_pre(&u[1]).unwrap()];

        let w_tv2 = inv2(&mul2(&pre[0].tv2, &pre[1].tv2));
        let invs = [mul2(&w_tv2, &pre[1].tv2), mul2(&w_tv2, &pre[0].tv2)];

        let mut flags = [0u8; 2];
        let mut ys = [Fp2 { c0: ZERO, c1: ZERO }; 2];
        let mut points = Vec::new();
        for i in 0..2 {
            let (flag, x, gx) = select_sswu_branch(&pre[i], &invs[i]);
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
        let lambda = mul2(&sub2(&points[1].y, &points[0].y), &inv2(&dx));
        let x3 = sub2(&sub2(&sq2(&lambda), &points[0].x), &points[1].x);
        let y3 = sub2(&mul2(&lambda, &sub2(&points[0].x, &x3)), &points[0].y);
        let (x_out, y_out) = iso_outputs(&Point2 { x: x3, y: y3 });

        let mut blob = Vec::with_capacity(W_TOTAL);
        blob.push(flags[0]);
        push_fp2_mont(&mut blob, &ys[0]);
        blob.push(flags[1]);
        push_fp2_mont(&mut blob, &ys[1]);
        #[cfg(not(feature = "wide-witness"))]
        push_fp2_mont(&mut blob, &w_tv2);
        #[cfg(feature = "wide-witness")]
        {
            push_fp2_mont(&mut blob, &invs[0]);
            push_fp2_mont(&mut blob, &invs[1]);
        }
        push_fp2_mont(&mut blob, &lambda);
        push_fp2_mont(&mut blob, &x_out);
        push_fp2_mont(&mut blob, &y_out);

        assert_eq!(blob.len(), W_TOTAL);
        blob
    }

    /// The 97-byte prefix (flags and the two root halves) plus the batch
    /// product, so each layout pays only for what it ships. The batch
    /// elements mirror compact_map_parts, which the verifier also uses.
    fn generate_compact_prefix(msg: &[u8], steer_flags: u8) -> (Vec<u8>, Fp) {
        let u = hash_to_field_g2(crate::dst::G2_RO, msg);
        let pre = [sswu_pre(&u[0]).unwrap(), sswu_pre(&u[1]).unwrap()];

        let mut flags = 0u8;
        let mut blob = Vec::with_capacity(C_TOTAL);
        blob.push(0);
        let mut prod = R;
        for i in 0..2 {
            let (mut flag, _, gx) = select_sswu_branch(&pre[i], &inv2(&pre[i].tv2));
            let y = sqrt2(&gx);
            assert!(
                !is_zero(&y.c0),
                "measure-zero input: the root's real part is zero, \
                 which the compact blob cannot encode"
            );
            // test-only steering: lie about the branch while keeping the
            // batch self-consistent with the lie
            flag ^= (steer_flags >> i) & 1;
            flags |= flag << i;
            blob.extend_from_slice(&limbs_to_be(&y.c0));

            // the product is order-independent, so just fold the four
            // per-map batch elements in as they appear
            let parts = compact_map_parts(&pre[i], flag);
            prod = mont_mul(&prod, &add_mod(&y.c0, &y.c0));
            prod = mont_mul(&prod, &norm_fp2(&pre[i].tv2));
            prod = mont_mul(&prod, &norm_fp2(&parts.xden_h));
            prod = mont_mul(&prod, &norm_fp2(&parts.yden_h));
        }
        blob[0] = flags;
        (blob, prod)
    }

    /// The 145-byte compact blob: flags, the two roots' real parts, one
    /// batched inverse.
    pub fn generate_compact(msg: &[u8]) -> Vec<u8> {
        generate_compact_steered(msg, 0)
    }

    /// The 97-byte witness-free-inverse blob: the compact blob minus its
    /// batched-inverse tail (the program recomputes that by binary xgcd).
    pub fn generate_compact_xgcd(msg: &[u8]) -> Vec<u8> {
        generate_compact_prefix(msg, 0).0
    }

    /// The other square root for the first map: an equally valid witness
    /// that must not steer. Negating c0 negates one batch factor, so the
    /// trailing inverse (when the layout carries it) negates with it.
    pub fn flip_first_root(blob: &[u8]) -> Vec<u8> {
        assert!(blob.len() == C_TOTAL || blob.len() == C_XGCD_TOTAL);
        let mut out = blob.to_vec();
        for range in [1..49, C_XGCD_TOTAL..C_TOTAL] {
            if range.end > blob.len() {
                break;
            }
            let flipped = neg_mod(&be_to_limbs(blob[range.clone()].try_into().unwrap()));
            out[range].copy_from_slice(&limbs_to_be(&flipped));
        }
        out
    }

    /// Branch-lie probe: the batch witness stays consistent with the lie,
    /// but no c0 can satisfy the real-part check against a non-square gx,
    /// so the program must abort.
    #[doc(hidden)]
    pub fn generate_compact_steered(msg: &[u8], steer_flags: u8) -> Vec<u8> {
        let (mut blob, prod) = generate_compact_prefix(msg, steer_flags);
        blob.extend_from_slice(&limbs_to_be(&inv_fp(&prod)));
        blob
    }

    /// Branch-lie probe, 97-byte layout.
    #[doc(hidden)]
    pub fn generate_compact_xgcd_steered(msg: &[u8], steer_flags: u8) -> Vec<u8> {
        generate_compact_prefix(msg, steer_flags).0
    }
}
