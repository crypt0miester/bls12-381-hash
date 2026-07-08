//! Witnessed Shallue-van de Woestijne hash to G2: maps field elements
//! directly onto E2, skipping the 3-isogeny of the SSWU pipeline. Map
//! construction from Wahby-Boneh eprint 2019/403 section 3 with u0 = -1
//! (their exception-free choice for E2). Same witness discipline and
//! branch-proof structure as g1_svdw, lifted to Fp2: non-squareness of
//! f(x_i) is proven with s^2 = (1+i) f(x_i), reusing the 1+i non-residue
//! already underpinning the SSWU branch trick.
//!
//! This is NOT the standard ciphersuite: output differs from
//! BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_ (blst, EIP-2537), and the
//! DST is changed accordingly.

use solana_program::program_error::ProgramError;

use crate::g1_msig::{add_mod, mont_mul, sub_mod};
use crate::g2_consts::{SVDW2_B, SVDW2_C1, SVDW2_F_U0, SVDW2_INV_3U0SQ, SVDW2_SQRT_M3};
use crate::g2_msig::{
    add2, add_prime_witnessed, check_inverse2, clear_cofactor, fold, fp2, from_mont2,
    g2_validate, is_zero2, mul2, neg2, point_bytes, sgn0_fp2, sq2, sub2, to_mont2, wit96,
    Elem2, Fp2, Point2, ONE2, ZERO,
};

const DST_SVDW2: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SVDW_RO_POP_";

/// Witness bytes for a map claiming branch j: flag, the shared inverse
/// alpha = 1/(t^2 (t^2 + f(u0))), j-1 non-squareness proofs, one sqrt.
fn map_witness_len(j: u8) -> usize {
    97 + 96 * j as usize
}

fn expand_message_xmd(msg: &[u8]) -> [[u8; 32]; 8] {
    use solana_program::hash::hashv;

    let z_pad = [0u8; 64];
    let l_i_b = [1u8, 0];
    let dst_len = [DST_SVDW2.len() as u8];

    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], DST_SVDW2, &dst_len]).to_bytes();

    let mut blocks = [[0u8; 32]; 8];
    blocks[0] = hashv(&[&b0, &[1u8], DST_SVDW2, &dst_len]).to_bytes();
    for i in 1..8 {
        let mut x = [0u8; 32];
        for j in 0..32 {
            x[j] = b0[j] ^ blocks[i - 1][j];
        }
        blocks[i] = hashv(&[&x, &[i as u8 + 1], DST_SVDW2, &dst_len]).to_bytes();
    }
    blocks
}

fn hash_to_field(msg: &[u8]) -> [Elem2; 2] {
    let blocks = expand_message_xmd(msg);
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

/// Multiply by a base-field constant: (a + bi) k = ak + bk i.
fn mul2_real(a: &Fp2, k: &[u64; 6]) -> Fp2 {
    Fp2 {
        c0: mont_mul(&a.c0, k),
        c1: mont_mul(&a.c1, k),
    }
}

/// Multiply by 1+i: (a + bi)(1 + i) = (a - b) + (a + b) i. Domain-free.
fn mul2_one_plus_i(a: &Fp2) -> Fp2 {
    Fp2 {
        c0: sub_mod(&a.c0, &a.c1),
        c1: add_mod(&a.c0, &a.c1),
    }
}

/// f(x) = x^3 + 4(1+i) on E2 (a = 0), Montgomery in and out.
pub(crate) fn f_at(x: &Fp2) -> Fp2 {
    add2(&mul2(&sq2(x), x), &fp2(&SVDW2_B))
}

/// Verify that fx is a non-residue: witness s with s^2 = (1+i) fx. Sound
/// because 1+i is a non-residue (its norm 2 is a non-square mod p), so
/// (1+i) fx square implies fx non-square. fx = 0 is rejected: E2 has odd
/// order, so f has no roots and the case is unreachable for honest
/// inputs, but a zero witness would otherwise steer the branch.
fn check_nonsquare(fx: &Fp2, wit: &[u8]) -> Result<(), ProgramError> {
    if is_zero2(fx) {
        return Err(ProgramError::InvalidInstructionData);
    }
    let s_m = wit96(wit)?;
    if sq2(&s_m) != mul2_one_plus_i(fx) {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(())
}

/// SvdW map on input t, branch and expensive results supplied as witness.
///
/// x-candidates per eprint 2019/403 section 3, all reachable from the one
/// inverse alpha = 1/(t^2 (t^2 + f(u0))):
///   x1 = (sqrt(-3) - u0)/2 - alpha t^4 sqrt(-3)
///   x2 = -u0 - x1
///   x3 = u0 - alpha (t^2 + f(u0))^3 / 3
fn map_to_curve_witnessed(t: &Elem2, wit: &[u8]) -> Result<Point2, ProgramError> {
    let j = wit[0];
    let w_alpha = wit96(&wit[1..97])?;

    let t2 = sq2(&t.mont);
    let t4 = sq2(&t2);
    let f_t = add2(&t2, &fp2(&SVDW2_F_U0));
    let q = mul2(&t2, &f_t);
    // q == 0 needs t = 0 or t^2 = -f(u0); probability ~2^-760, reject
    if is_zero2(&q) {
        return Err(ProgramError::InvalidInstructionData);
    }
    let alpha = check_inverse2(&q, &w_alpha)?;

    let neg_u0 = ONE2;
    let u0 = neg2(&neg_u0);
    let x1 = sub2(
        &Fp2 { c0: SVDW2_C1, c1: ZERO },
        &mul2_real(&mul2(&t4, &alpha), &SVDW2_SQRT_M3),
    );

    let mut off = 97;
    let x = match j {
        1 => x1,
        2 => {
            check_nonsquare(&f_at(&x1), &wit[off..off + 96])?;
            off += 96;
            sub2(&neg_u0, &x1)
        }
        3 => {
            let x2 = sub2(&neg_u0, &x1);
            check_nonsquare(&f_at(&x1), &wit[off..off + 96])?;
            off += 96;
            check_nonsquare(&f_at(&x2), &wit[off..off + 96])?;
            off += 96;
            let ft3 = mul2(&sq2(&f_t), &f_t);
            sub2(&u0, &mul2_real(&mul2(&ft3, &alpha), &SVDW2_INV_3U0SQ))
        }
        _ => return Err(ProgramError::InvalidInstructionData),
    };

    let gx = f_at(&x);
    let y_w = wit96(&wit[off..off + 96])?;
    let yw_m = to_mont2(&y_w);
    if sq2(&yw_m) != gx {
        return Err(ProgramError::InvalidInstructionData);
    }

    // sgn0 correction: sign of y must match sign of t, which also pins
    // the output against the witness supplying the other root
    let mut y_canonical = y_w;
    if sgn0_fp2(&y_canonical) != sgn0_fp2(&t.canonical) {
        y_canonical = neg2(&y_canonical);
    }

    Ok(Point2 { x, y: to_mont2(&y_canonical) })
}

/// Payload: [map0 witness][map1 witness][dx inverse: 96][message].
pub fn run_witnessed(payload: &[u8]) -> Result<Vec<u8>, ProgramError> {
    let split_map = |data: &[u8]| -> Result<usize, ProgramError> {
        let j = *data.first().ok_or(ProgramError::InvalidInstructionData)?;
        if !(1..=3).contains(&j) {
            return Err(ProgramError::InvalidInstructionData);
        }
        let len = map_witness_len(j);
        if data.len() < len {
            return Err(ProgramError::InvalidInstructionData);
        }
        Ok(len)
    };

    let len0 = split_map(payload)?;
    let (wit0, rest) = payload.split_at(len0);
    let len1 = split_map(rest)?;
    let (wit1, rest) = rest.split_at(len1);
    if rest.len() < 96 {
        return Err(ProgramError::InvalidInstructionData);
    }
    let (dx_bytes, msg) = rest.split_at(96);
    let w_dx = wit96(dx_bytes)?;

    let u = hash_to_field(msg);
    let p0 = map_to_curve_witnessed(&u[0], wit0)?;
    let p1 = map_to_curve_witnessed(&u[1], wit1)?;
    let sum = add_prime_witnessed(&p0, &p1, &w_dx)?;

    let uncleared = point_bytes(&from_mont2(&sum.x), &from_mont2(&sum.y));
    let cleared = clear_cofactor(&uncleared)?;
    g2_validate(&cleared)?;
    Ok(cleared.to_vec())
}

/// Host-side witness generation and reference evaluation, mirroring the
/// on-chain pipeline with the expensive results computed via
/// square-and-multiply.
#[cfg(not(target_os = "solana"))]
pub mod witness {
    use super::*;
    use crate::g2_msig::witness::{inv2, is_square2, push_fp2, push_fp2_mont, sqrt2};

    /// Montgomery-form mapped point plus its serialized witness blob.
    fn map_one(t: &Elem2) -> (Point2, Vec<u8>) {
        let t2 = sq2(&t.mont);
        let t4 = sq2(&t2);
        let f_t = add2(&t2, &fp2(&SVDW2_F_U0));
        let q = mul2(&t2, &f_t);
        assert!(!is_zero2(&q));

        let alpha = inv2(&q);

        let neg_u0 = ONE2;
        let u0 = neg2(&neg_u0);
        let x1 = sub2(
            &Fp2 { c0: SVDW2_C1, c1: ZERO },
            &mul2_real(&mul2(&t4, &alpha), &SVDW2_SQRT_M3),
        );
        let x2 = sub2(&neg_u0, &x1);
        let ft3 = mul2(&sq2(&f_t), &f_t);
        let x3 = sub2(&u0, &mul2_real(&mul2(&ft3, &alpha), &SVDW2_INV_3U0SQ));

        let mut blob = Vec::new();
        let mut chosen = None;
        for (j, x) in [(1u8, x1), (2, x2), (3, x3)] {
            let fx = f_at(&x);
            if is_square2(&fx) {
                chosen = Some((j, x, fx));
                break;
            }
            // non-squareness proof: sqrt of (1+i) fx
            let s = sqrt2(&mul2_one_plus_i(&fx));
            push_fp2_mont(&mut blob, &s);
        }
        let (j, x, fx) = chosen.expect("SvdW: one of three candidates is square");

        let y = sqrt2(&fx);
        push_fp2(&mut blob, &y);

        let mut header = vec![j];
        push_fp2_mont(&mut header, &alpha);
        header.extend_from_slice(&blob);
        assert_eq!(header.len(), map_witness_len(j));

        let mut y_canonical = from_mont2(&y);
        if sgn0_fp2(&y_canonical) != sgn0_fp2(&t.canonical) {
            y_canonical = neg2(&y_canonical);
        }
        (Point2 { x, y: to_mont2(&y_canonical) }, header)
    }

    fn mapped_points(msg: &[u8]) -> ([Point2; 2], Vec<u8>) {
        let u = hash_to_field(msg);
        let (p0, blob0) = map_one(&u[0]);
        let (p1, blob1) = map_one(&u[1]);
        let mut blob = blob0;
        blob.extend_from_slice(&blob1);
        ([p0, p1], blob)
    }

    pub fn generate(msg: &[u8]) -> Vec<u8> {
        let ([p0, p1], mut blob) = mapped_points(msg);
        let dx = sub2(&p1.x, &p0.x);
        push_fp2_mont(&mut blob, &inv2(&dx));
        blob
    }

    /// The sum point before cofactor clearing, uncompressed affine bytes.
    /// Callers apply the effective cofactor to get the expected output.
    pub fn reference_preclear(msg: &[u8]) -> [u8; 192] {
        let ([p0, p1], _) = mapped_points(msg);
        let dx = sub2(&p1.x, &p0.x);
        let lambda = mul2(&sub2(&p1.y, &p0.y), &inv2(&dx));
        let x3 = sub2(&sub2(&sq2(&lambda), &p0.x), &p1.x);
        let y3 = sub2(&mul2(&lambda, &sub2(&p0.x, &x3)), &p0.y);
        point_bytes(&from_mont2(&x3), &from_mont2(&y3))
    }

    /// Replace map0's sqrt witness with the other root. A valid witness
    /// either way; the on-chain sign rule must produce the same point.
    pub fn flip_first_sqrt(blob: &[u8]) -> Vec<u8> {
        let j = blob[0];
        let y_off = map_witness_len(j) - 96;
        let y = wit96(&blob[y_off..y_off + 96]).unwrap();
        let mut out = blob[..y_off].to_vec();
        push_fp2(&mut out, &to_mont2(&neg2(&y)));
        out.extend_from_slice(&blob[y_off + 96..]);
        out
    }
}
