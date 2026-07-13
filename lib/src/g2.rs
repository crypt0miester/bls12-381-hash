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
    add_mod, be_to_limbs, limbs_to_be, mont_mul, neg_mod, sub_mod, sys,
    to_mont, Fp,
};
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

/// Adapted iso-3 evaluation: cubics as (y + g)(w + q1) + q0, the degree-2
/// denominator via its tiny coefficients (12 - 12i, -72i) as adds. Five
/// Fp2 multiplications plus one squaring against eleven for Horner.
pub(crate) fn iso3_adapted(x: &Fp2) -> (Fp2, Fp2, Fp2, Fp2) {
    fn m12(a: &Fp) -> Fp {
        let a2 = add_mod(a, a);
        let a4 = add_mod(&a2, &a2);
        add_mod(&add_mod(&a4, &a4), &a4)
    }
    let w = sq2(x);

    // sums ride unreduced: every one is a mul2 operand or checked via mul2
    let cubic = |k: &[[[u64; 6]; 2]]| {
        add2_unreduced(
            &mul2(&add2_unreduced(x, &fp2(&k[0])), &add2_unreduced(&w, &fp2(&k[1]))),
            &fp2(&k[2]),
        )
    };
    let x_num = mul2(&cubic(&ISO3A_XNUM), &fp2(&ISO3A_XNUM[3]));
    let y_num = mul2(&cubic(&ISO3A_YNUM), &fp2(&ISO3A_YNUM[3]));
    let y_den = cubic(&ISO3A_YDEN);
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

/// Host-side witness generation mirroring the on-chain pipeline.
#[cfg(not(target_os = "solana"))]
pub mod witness {
    use super::*;
    use crate::consts_g1::R;
    use crate::consts_g2::{ISO3_XDEN, ISO3_XNUM, ISO3_YDEN, ISO3_YNUM};
    use crate::fp::{add_carryless, exp_inverse, exp_legendre, exp_sqrt, is_zero, shr1};
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
        push_fp2_mont(&mut out, &neg2(&y));
        out.extend_from_slice(&blob[97..]);
        out
    }

    pub fn generate_nu(msg: &[u8]) -> Vec<u8> {
        let u = hash_to_field_nu(crate::dst::G2_NU, msg);
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
}
