//! Witnessed Shallue-van de Woestijne hash to G1: maps field elements
//! directly onto E1, so the 11-isogeny evaluation that dominates the SSWU
//! pipeline disappears entirely. Map construction from Wahby-Boneh eprint
//! 2019/403 section 3 with u0 = -3 (their exception-free choice for E1).
//!
//! This is NOT the standard ciphersuite: output differs from
//! BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_ (blst, EIP-2537), and the
//! DST is changed accordingly. Interop with existing BLS stacks is the
//! price of skipping the isogeny.

use solana_program_error::ProgramError;
use alloc::vec::Vec;

use crate::consts_g2::C256_MONT;
use crate::consts_g1::{
    SVDW_B, SVDW_C1, SVDW_F_U0, SVDW_INV_3U0SQ, SVDW_NEG_U0, SVDW_SQRT_M27, SVDW_U0,
};
use crate::fp::*;
use crate::g1::*;

const DST_SVDW: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SVDW_RO_POP_";

/// Witness bytes for a map claiming branch j: flag, the shared inverse
/// alpha = 1/(t^2 (t^2 + f(u0))), j-1 non-squareness proofs, one sqrt.
fn map_witness_len(j: u8) -> usize {
    49 + 48 * j as usize
}

struct FieldElem {
    canonical: Fp,
    mont: Fp,
}

/// Affine point on E1, Montgomery form coordinates.
struct PointE1 {
    x: Fp,
    y: Fp,
}

fn expand_message_xmd(msg: &[u8]) -> [[u8; 32]; 4] {
    use solana_sha256_hasher::hashv;

    let z_pad = [0u8; 64];
    let l_i_b = [0u8, 128];
    let dst_len = [DST_SVDW.len() as u8];

    let b0 = hashv(&[&z_pad, msg, &l_i_b, &[0u8], DST_SVDW, &dst_len]).to_bytes();

    let mut blocks = [[0u8; 32]; 4];
    blocks[0] = hashv(&[&b0, &[1u8], DST_SVDW, &dst_len]).to_bytes();
    for i in 1..4 {
        let mut x = [0u8; 32];
        for j in 0..32 {
            x[j] = b0[j] ^ blocks[i - 1][j];
        }
        blocks[i] = hashv(&[&x, &[i as u8 + 1], DST_SVDW, &dst_len]).to_bytes();
    }
    blocks
}

/// hash_to_field without modexp: split the 64-byte value at bit 256 and
/// fold with a precomputed 2^256 mod p.
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

/// f(x) = x^3 + 4 on E1 (a = 0), Montgomery in and out.
fn f_at(x: &Fp) -> Fp {
    add_mod(&mont_mul(&mont_sqr(x), x), &SVDW_B)
}

/// Verify that fx is a non-residue: witness s with s^2 = -fx. 
fn check_nonsquare(fx: &Fp, wit: &[u8]) -> Result<(), ProgramError> {
    if is_zero(fx) {
        return Err(ProgramError::InvalidInstructionData);
    }
    let s_m = wit48(wit)?;
    if mont_sqr(&s_m) != neg_mod(fx) {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(())
}

/// SvdW map on input t, branch and expensive results supplied as witness.
fn map_to_curve_witnessed(t: &FieldElem, wit: &[u8]) -> Result<PointE1, ProgramError> {
    let j = wit[0];
    let w_alpha = wit48(&wit[1..49])?;

    let t2 = mont_sqr(&t.mont);
    let t4 = mont_sqr(&t2);
    let f_t = add_mod(&t2, &SVDW_F_U0);
    let q = mont_mul(&t2, &f_t);
    // q == 0 needs t = 0 or t^2 = 23; probability ~2^-380, reject
    if is_zero(&q) {
        return Err(ProgramError::InvalidInstructionData);
    }
    let alpha_m = check_inverse(&q, &w_alpha)?;

    let x1 = sub_mod(
        &SVDW_C1,
        &mont_mul(&mont_mul(&t4, &alpha_m), &SVDW_SQRT_M27),
    );

    let mut off = 49;
    let x = match j {
        1 => x1,
        2 => {
            check_nonsquare(&f_at(&x1), &wit[off..off + 48])?;
            off += 48;
            sub_mod(&SVDW_NEG_U0, &x1)
        }
        3 => {
            let x2 = sub_mod(&SVDW_NEG_U0, &x1);
            check_nonsquare(&f_at(&x1), &wit[off..off + 48])?;
            off += 48;
            check_nonsquare(&f_at(&x2), &wit[off..off + 48])?;
            off += 48;
            let ft3 = mont_mul(&mont_sqr(&f_t), &f_t);
            sub_mod(
                &SVDW_U0,
                &mont_mul(&mont_mul(&ft3, &alpha_m), &SVDW_INV_3U0SQ),
            )
        }
        _ => return Err(ProgramError::InvalidInstructionData),
    };

    let gx = f_at(&x);
    let y_w = wit48(&wit[off..off + 48])?;
    let yw_m = to_mont(&y_w);
    if mont_sqr(&yw_m) != gx {
        return Err(ProgramError::InvalidInstructionData);
    }

    // sgn0 correction: parity of y must match parity of t, which also
    // pins the output against the witness supplying the other root
    let mut y_canonical = y_w;
    if (y_canonical[0] & 1) != (t.canonical[0] & 1) {
        y_canonical = neg_mod(&y_canonical);
    }

    Ok(PointE1 { x, y: to_mont(&y_canonical) })
}

/// Affine chord addition on E1 (coefficient-free; errors on x1 == x2).
fn add_e1_witnessed(p: &PointE1, q: &PointE1, w_dx: &Fp) -> Result<PointE1, ProgramError> {
    if p.x == q.x {
        return Err(ProgramError::InvalidInstructionData);
    }
    let dx = sub_mod(&q.x, &p.x);
    let inv_m = check_inverse(&dx, w_dx)?;
    let lambda = mont_mul(&sub_mod(&q.y, &p.y), &inv_m);
    let x3 = sub_mod(&sub_mod(&mont_sqr(&lambda), &p.x), &q.x);
    let y3 = sub_mod(&mont_mul(&lambda, &sub_mod(&p.x, &x3)), &p.y);
    Ok(PointE1 { x: x3, y: y3 })
}

/// Payload: [map0 witness][map1 witness][dx inverse: 48][message].
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
    if rest.len() < 48 {
        return Err(ProgramError::InvalidInstructionData);
    }
    let (dx_bytes, msg) = rest.split_at(48);
    let w_dx = wit48(dx_bytes)?;

    let u = hash_to_field_folded(msg);
    let p0 = map_to_curve_witnessed(&u[0], wit0)?;
    let p1 = map_to_curve_witnessed(&u[1], wit1)?;
    let sum = add_e1_witnessed(&p0, &p1, &w_dx)?;

    let x = limbs_to_be(&from_mont(&sum.x));
    let y = limbs_to_be(&from_mont(&sum.y));
    let cleared = clear_cofactor(&point_bytes(&x, &y))?;
    validate(&cleared)?;
    Ok(cleared.to_vec())
}

/// Host-side witness generation and reference evaluation, mirroring the
/// on-chain pipeline with the expensive results computed via
/// square-and-multiply.
#[cfg(not(target_os = "solana"))]
pub mod witness {
    use super::*;
    use crate::consts_g1::R;
    
    

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

    /// Montgomery-form inverse witness bytes.
    fn inverse_witness(v_m: &Fp) -> [u8; 48] {
        limbs_to_be(&pow_mont(v_m, &exp_inverse()))
    }

    fn is_square(v_m: &Fp) -> bool {
        is_zero(v_m) || pow_mont(v_m, &exp_legendre()) == R
    }

    /// Montgomery-form mapped point plus its serialized witness blob.
    fn map_one(t: &FieldElem) -> (PointE1, Vec<u8>) {
        let t2 = mont_sqr(&t.mont);
        let t4 = mont_sqr(&t2);
        let f_t = add_mod(&t2, &SVDW_F_U0);
        let q = mont_mul(&t2, &f_t);
        assert!(!is_zero(&q));

        let w_alpha = inverse_witness(&q);
        let alpha_m = be_to_limbs(&w_alpha);

        let x1 = sub_mod(
            &SVDW_C1,
            &mont_mul(&mont_mul(&t4, &alpha_m), &SVDW_SQRT_M27),
        );
        let x2 = sub_mod(&SVDW_NEG_U0, &x1);
        let ft3 = mont_mul(&mont_sqr(&f_t), &f_t);
        let x3 = sub_mod(
            &SVDW_U0,
            &mont_mul(&mont_mul(&ft3, &alpha_m), &SVDW_INV_3U0SQ),
        );

        let mut blob = Vec::new();
        let mut chosen = None;
        for (j, x) in [(1u8, x1), (2, x2), (3, x3)] {
            let fx = f_at(&x);
            if is_square(&fx) {
                chosen = Some((j, x, fx));
                break;
            }
            // non-squareness proof: sqrt of -fx
            let s_m = pow_mont(&neg_mod(&fx), &exp_sqrt());
            assert_eq!(mont_sqr(&s_m), neg_mod(&fx));
            blob.extend_from_slice(&limbs_to_be(&s_m));
        }
        let (j, x, fx) = chosen.expect("SvdW: one of three candidates is square");

        let y_m = pow_mont(&fx, &exp_sqrt());
        assert_eq!(mont_mul(&y_m, &y_m), fx);
        let y_c = from_mont(&y_m);
        blob.extend_from_slice(&limbs_to_be(&y_c));

        let mut header = vec![j];
        header.extend_from_slice(&w_alpha);
        header.extend_from_slice(&blob);
        assert_eq!(header.len(), map_witness_len(j));

        let mut y_final = y_c;
        if (y_final[0] & 1) != (t.canonical[0] & 1) {
            y_final = neg_mod(&y_final);
        }
        (PointE1 { x, y: to_mont(&y_final) }, header)
    }

    fn mapped_points(msg: &[u8]) -> ([PointE1; 2], Vec<u8>) {
        let u = hash_to_field_folded(msg);
        let (p0, blob0) = map_one(&u[0]);
        let (p1, blob1) = map_one(&u[1]);
        let mut blob = blob0;
        blob.extend_from_slice(&blob1);
        ([p0, p1], blob)
    }

    pub fn generate(msg: &[u8]) -> Vec<u8> {
        let ([p0, p1], mut blob) = mapped_points(msg);
        let dx = sub_mod(&p1.x, &p0.x);
        blob.extend_from_slice(&inverse_witness(&dx));
        blob
    }

    /// The sum point before cofactor clearing, uncompressed affine bytes.
    /// Callers apply the effective cofactor to get the expected output.
    pub fn reference_preclear(msg: &[u8]) -> [u8; 96] {
        let ([p0, p1], _) = mapped_points(msg);
        let dx = sub_mod(&p1.x, &p0.x);
        let inv_m = be_to_limbs(&inverse_witness(&dx));
        let lambda = mont_mul(&sub_mod(&p1.y, &p0.y), &inv_m);
        let x3 = sub_mod(&sub_mod(&mont_sqr(&lambda), &p0.x), &p1.x);
        let y3 = sub_mod(&mont_mul(&lambda, &sub_mod(&p0.x, &x3)), &p0.y);
        point_bytes(
            &limbs_to_be(&from_mont(&x3)),
            &limbs_to_be(&from_mont(&y3)),
        )
    }

    /// Replace map0's sqrt witness with the other root. A valid witness
    /// either way; the on-chain parity rule must produce the same point.
    pub fn flip_first_sqrt(blob: &[u8]) -> Vec<u8> {
        let j = blob[0];
        let y_off = map_witness_len(j) - 48;
        let y = be_to_limbs(blob[y_off..y_off + 48].try_into().unwrap());
        let mut out = blob.to_vec();
        out[y_off..y_off + 48].copy_from_slice(&limbs_to_be(&neg_mod(&y)));
        out
    }
}
