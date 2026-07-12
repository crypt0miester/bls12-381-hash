#![allow(unexpected_cfgs)]

use core::hint::black_box;

use bls381_hash::dst::{G1_NU, G1_RO, G2_NU, G2_RO};
use bls381_hash::{
    encode_to_g1, encode_to_g2, hash_to_g1, hash_to_g1_modexp, hash_to_g2,
};

use bls12_381::{
    hash_to_curve::{ExpandMsgXmd, HashToCurve, HashToField, MapToCurve},
    G1Affine, G1Projective, G2Affine, G2Projective,
};
use sha2::Sha256;
use solana_program::{
    account_info::AccountInfo,
    entrypoint,
    entrypoint::ProgramResult,
    program::set_return_data,
    program_error::ProgramError,
    pubkey::Pubkey,
};


const BLS12_381_PAIRING_BE: u64 = 4 | 0x80;
const BLS12_381_G1_BE: u64 = 5 | 0x80;
const BLS12_381_G2_BE: u64 = 6 | 0x80;

const OP_ADD: u64 = 0;
const OP_SUB: u64 = 1;
const OP_MUL: u64 = 2;

const G1_POINT: usize = 96;
const G1_COMPRESSED: usize = 48;
const G2_POINT: usize = 192;
const G2_COMPRESSED: usize = 96;
const SCALAR: usize = 32;
const GT: usize = 576;

#[cfg(target_os = "solana")]
mod sys {
    use solana_define_syscall::define_syscall;

    define_syscall!(fn sol_curve_validate_point(curve_id: u64, point_addr: *const u8, result: *mut u8) -> u64);
    define_syscall!(fn sol_curve_group_op(curve_id: u64, group_op: u64, left_input_addr: *const u8, right_input_addr: *const u8, result_point_addr: *mut u8) -> u64);
    define_syscall!(fn sol_curve_pairing_map(curve_id: u64, num_pairs: u64, g1_points: *const u8, g2_points: *const u8, result: *mut u8) -> u64);
    define_syscall!(fn sol_curve_decompress(curve_id: u64, point: *const u8, result: *mut u8) -> u64);
    define_syscall!(fn sol_big_mod_exp(params: *const u8, result: *mut u8) -> u64);
    define_syscall!(fn sol_memcmp_(s1: *const u8, s2: *const u8, n: u64, result: *mut i32) -> u64);
}

#[cfg(not(target_os = "solana"))]
#[allow(clippy::missing_safety_doc)]
mod sys {
    pub unsafe fn sol_curve_validate_point(_: u64, _: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_curve_group_op(_: u64, _: u64, _: *const u8, _: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_curve_pairing_map(_: u64, _: u64, _: *const u8, _: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_curve_decompress(_: u64, _: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_big_mod_exp(_: *const u8, _: *mut u8) -> u64 {
        unimplemented!()
    }
    pub unsafe fn sol_memcmp_(_: *const u8, _: *const u8, _: u64, _: *mut i32) -> u64 {
        unimplemented!()
    }
}

/// The GT identity in the pairing syscall's encoding: big-endian Fp12
/// coefficients, so the constant one sits in the trailing byte.
const GT_ONE: [u8; GT] = {
    let mut one = [0u8; GT];
    one[GT - 1] = 1;
    one
};

/// The negated G1 generator, uncompressed: same x, y = p - y_gen. Baked in
/// so the two-pair product check needs no generator in the payload.
const NEG_G1_GEN: [u8; 96] = [
    0x17, 0xf1, 0xd3, 0xa7, 0x31, 0x97, 0xd7, 0x94, 0x26, 0x95, 0x63, 0x8c,
    0x4f, 0xa9, 0xac, 0x0f, 0xc3, 0x68, 0x8c, 0x4f, 0x97, 0x74, 0xb9, 0x05,
    0xa1, 0x4e, 0x3a, 0x3f, 0x17, 0x1b, 0xac, 0x58, 0x6c, 0x55, 0xe8, 0x3f,
    0xf9, 0x7a, 0x1a, 0xef, 0xfb, 0x3a, 0xf0, 0x0a, 0xdb, 0x22, 0xc6, 0xbb,
    0x11, 0x4d, 0x1d, 0x68, 0x55, 0xd5, 0x45, 0xa8, 0xaa, 0x7d, 0x76, 0xc8,
    0xcf, 0x2e, 0x21, 0xf2, 0x67, 0x81, 0x6a, 0xef, 0x1d, 0xb5, 0x07, 0xc9,
    0x66, 0x55, 0xb9, 0xd5, 0xca, 0xac, 0x42, 0x36, 0x4e, 0x6f, 0x38, 0xba,
    0x0e, 0xcb, 0x75, 0x1b, 0xad, 0x54, 0xdc, 0xd6, 0xb9, 0x39, 0xc2, 0xca,
];

#[repr(C)]
struct BigModExpParams {
    base: *const u8,
    base_len: u64,
    exponent: *const u8,
    exponent_len: u64,
    modulus: *const u8,
    modulus_len: u64,
}

entrypoint!(process_instruction);

fn process_instruction(
    _program_id: &Pubkey,
    _accounts: &[AccountInfo],
    data: &[u8],
) -> ProgramResult {
    let (&tag, payload) = data
        .split_first()
        .ok_or(ProgramError::InvalidInstructionData)?;

    match tag {
        // In-program hash-to-curve pipeline, cumulative prefixes so the host
        // can subtract successive results to isolate each phase.
        0 => {
            let p = <G2Projective as HashToCurve<ExpandMsgXmd<Sha256>>>::hash_to_curve(
                payload, G2_RO,
            );
            set_return_data(&G2Affine::from(p).to_compressed());
        }
        1 => {
            let p = <G1Projective as HashToCurve<ExpandMsgXmd<Sha256>>>::hash_to_curve(
                payload, G1_RO,
            );
            set_return_data(&G1Affine::from(p).to_compressed());
        }
        2 => {
            let u = hash_to_field_g2(payload);
            black_box(&u);
        }
        3 => {
            let u = hash_to_field_g2(payload);
            let p = G2Projective::map_to_curve(&u[0]) + G2Projective::map_to_curve(&u[1]);
            black_box(&p);
        }
        4 => {
            let u = hash_to_field_g2(payload);
            let p = (G2Projective::map_to_curve(&u[0]) + G2Projective::map_to_curve(&u[1]))
                .clear_h();
            black_box(&p);
        }

        // Raw syscall probes; operands come from instruction data, results go
        // back through return data with the syscall status in the first byte.
        10 => {
            expect_len(payload, G2_POINT)?;
            let mut out = 0u8;
            let rc = unsafe {
                sys::sol_curve_validate_point(BLS12_381_G2_BE, payload.as_ptr(), &mut out)
            };
            set_return_data(&[rc as u8]);
        }
        11 => {
            expect_len(payload, 2 * G2_POINT)?;
            let mut out = [0u8; G2_POINT];
            let rc = unsafe {
                sys::sol_curve_group_op(
                    BLS12_381_G2_BE,
                    OP_ADD,
                    payload.as_ptr(),
                    payload[G2_POINT..].as_ptr(),
                    out.as_mut_ptr(),
                )
            };
            return_status(rc, &out);
        }
        12 => {
            expect_len(payload, SCALAR + G2_POINT)?;
            let mut out = [0u8; G2_POINT];
            let rc = unsafe {
                sys::sol_curve_group_op(
                    BLS12_381_G2_BE,
                    OP_MUL,
                    payload.as_ptr(),
                    payload[SCALAR..].as_ptr(),
                    out.as_mut_ptr(),
                )
            };
            return_status(rc, &out);
        }
        13 => {
            expect_len(payload, G2_COMPRESSED)?;
            let mut out = [0u8; G2_POINT];
            let rc = unsafe {
                sys::sol_curve_decompress(BLS12_381_G2_BE, payload.as_ptr(), out.as_mut_ptr())
            };
            return_status(rc, &out);
        }
        14 => {
            let pairs = payload.len() / (G1_POINT + G2_POINT);
            expect_len(payload, pairs * (G1_POINT + G2_POINT))?;
            let (g1s, g2s) = payload.split_at(pairs * G1_POINT);
            let mut out = [0u8; GT];
            let rc = unsafe {
                sys::sol_curve_pairing_map(
                    BLS12_381_PAIRING_BE,
                    pairs as u64,
                    g1s.as_ptr(),
                    g2s.as_ptr(),
                    out.as_mut_ptr(),
                )
            };
            return_status(rc, &out);
        }
        15 => {
            expect_len(payload, 2 * G1_POINT)?;
            let mut out = [0u8; G1_POINT];
            let rc = unsafe {
                sys::sol_curve_group_op(
                    BLS12_381_G1_BE,
                    OP_ADD,
                    payload.as_ptr(),
                    payload[G1_POINT..].as_ptr(),
                    out.as_mut_ptr(),
                )
            };
            return_status(rc, &out);
        }
        16 => {
            expect_len(payload, SCALAR + G1_POINT)?;
            let mut out = [0u8; G1_POINT];
            let rc = unsafe {
                sys::sol_curve_group_op(
                    BLS12_381_G1_BE,
                    OP_MUL,
                    payload.as_ptr(),
                    payload[SCALAR..].as_ptr(),
                    out.as_mut_ptr(),
                )
            };
            return_status(rc, &out);
        }
        17 => {
            expect_len(payload, G1_COMPRESSED)?;
            let mut out = [0u8; G1_POINT];
            let rc = unsafe {
                sys::sol_curve_decompress(BLS12_381_G1_BE, payload.as_ptr(), out.as_mut_ptr())
            };
            return_status(rc, &out);
        }
        18 => {
            expect_len(payload, G1_POINT)?;
            let mut out = 0u8;
            let rc = unsafe {
                sys::sol_curve_validate_point(BLS12_381_G1_BE, payload.as_ptr(), &mut out)
            };
            set_return_data(&[rc as u8]);
        }
        20 => {
            expect_len(payload, 3 * 48)?;
            let params = BigModExpParams {
                base: payload.as_ptr(),
                base_len: 48,
                exponent: payload[48..].as_ptr(),
                exponent_len: 48,
                modulus: payload[96..].as_ptr(),
                modulus_len: 48,
            };
            let mut out = [0u8; 48];
            let rc = unsafe {
                sys::sol_big_mod_exp(
                    &params as *const BigModExpParams as *const u8,
                    out.as_mut_ptr(),
                )
            };
            return_status(rc, &out);
        }
        // Syscall-assisted min-sig hash_to_G1, cumulative stage prefixes.
        30..=33 => {
            let out = hash_to_g1_modexp(G1_RO, tag - 30, payload)?;
            set_return_data(&out);
        }
        // End-to-end min-pk vote verify: witness hash_to_G2, subtract the
        // absentees from the stored committee aggregate, then check the
        // pairing equation as one two-pair product,
        // e(effective_pk, H) * e(-g1_gen, sig) == 1 in GT. The pairs share
        // the final exponentiation that two separate calls each pay (~13k
        // CU), the 576 byte identity compare rides the memcmp syscall, and
        // every syscall-filled buffer skips its zero-init.
        51 => {
            let absent = payload[0] as usize;
            let agg_end = 1 + G1_POINT;
            let sig_end = agg_end + G2_COMPRESSED;
            let abs_end = sig_end + 48 * absent;
            if payload.len() < abs_end {
                return Err(ProgramError::InvalidInstructionData);
            }

            let mut effective: [u8; G1_POINT] = payload[1..agg_end].try_into().unwrap();
            for i in 0..absent {
                let compressed = &payload[sig_end + i * 48..sig_end + (i + 1) * 48];
                let mut member = core::mem::MaybeUninit::<[u8; G1_POINT]>::uninit();
                let rc = unsafe {
                    sys::sol_curve_decompress(
                        BLS12_381_G1_BE,
                        compressed.as_ptr(),
                        member.as_mut_ptr() as *mut u8,
                    )
                };
                if rc != 0 {
                    return Err(ProgramError::InvalidInstructionData);
                }
                let mut out = core::mem::MaybeUninit::<[u8; G1_POINT]>::uninit();
                let rc = unsafe {
                    sys::sol_curve_group_op(
                        BLS12_381_G1_BE,
                        OP_SUB,
                        effective.as_ptr(),
                        member.as_ptr() as *const u8,
                        out.as_mut_ptr() as *mut u8,
                    )
                };
                if rc != 0 {
                    return Err(ProgramError::InvalidInstructionData);
                }
                // SAFETY: rc == 0 means the syscall wrote the whole point
                effective = unsafe { out.assume_init() };
            }

            let hash = hash_to_g2(G2_RO, &payload[abs_end..])?;

            // g1 and g2 sides of the two pairs, contiguous for the syscall;
            // the signature decompresses straight into its slot
            let mut g1s = [0u8; 2 * G1_POINT];
            g1s[..G1_POINT].copy_from_slice(&effective);
            g1s[G1_POINT..].copy_from_slice(&NEG_G1_GEN);

            let mut g2s = core::mem::MaybeUninit::<[u8; 2 * G2_POINT]>::uninit();
            let g2s_ptr = g2s.as_mut_ptr() as *mut u8;
            unsafe { core::ptr::copy_nonoverlapping(hash.as_ptr(), g2s_ptr, G2_POINT) };
            let rc = unsafe {
                sys::sol_curve_decompress(
                    BLS12_381_G2_BE,
                    payload[agg_end..sig_end].as_ptr(),
                    g2s_ptr.add(G2_POINT),
                )
            };
            if rc != 0 {
                return Err(ProgramError::InvalidInstructionData);
            }

            let mut gt = core::mem::MaybeUninit::<[u8; GT]>::uninit();
            let rc = unsafe {
                sys::sol_curve_pairing_map(
                    BLS12_381_PAIRING_BE,
                    2,
                    g1s.as_ptr(),
                    g2s.as_ptr() as *const u8,
                    gt.as_mut_ptr() as *mut u8,
                )
            };
            if rc != 0 {
                return Err(ProgramError::InvalidInstructionData);
            }

            let mut cmp = 0i32;
            unsafe {
                sys::sol_memcmp_(
                    gt.as_ptr() as *const u8,
                    GT_ONE.as_ptr(),
                    GT as u64,
                    &mut cmp,
                )
            };
            if cmp != 0 {
                return Err(ProgramError::InvalidInstructionData);
            }
            set_return_data(&[1]);
        }
        40 => {
            let out = hash_to_g1(G1_RO, payload)?;
            set_return_data(&out);
        }
        41 => {
            let out = hash_to_g2(G2_RO, payload)?;
            set_return_data(&out);
        }
        // Witnessed pipeline stage prefixes: payload is stage byte, blob, msg.
        46 => {
            let (&stage, rest) = payload
                .split_first()
                .ok_or(ProgramError::InvalidInstructionData)?;
            let out = bls381_hash::hash_to_g1_prefix(G1_RO, stage, rest)?;
            set_return_data(&out);
        }
        47 => {
            let (&stage, rest) = payload
                .split_first()
                .ok_or(ProgramError::InvalidInstructionData)?;
            let out = bls381_hash::hash_to_g2_prefix(G2_RO, stage, rest)?;
            set_return_data(&out);
        }
        // Witnessed encode_to_curve (RFC 9380 NU suites): single map.
        44 => {
            let out = encode_to_g1(G1_NU, payload)?;
            set_return_data(&out);
        }
        45 => {
            let out = encode_to_g2(G2_NU, payload)?;
            set_return_data(&out);
        }
        21 => {
            expect_len(payload, 8)?;
            let count = u64::from_le_bytes(payload.try_into().unwrap());
            let mut acc = 0u128;
            let mut x = 0x9e3779b97f4a7c15u64;
            for _ in 0..count {
                acc = acc.wrapping_add((x as u128).wrapping_mul(black_box(x) as u128));
                x = x.wrapping_add(0x6a09e667f3bcc909);
            }
            black_box(acc);
        }
        // Per function CU probe: payload is probe id then iteration count.
        24 => {
            expect_len(payload, 9)?;
            let id = payload[0];
            let count = u64::from_le_bytes(payload[1..9].try_into().unwrap());
            black_box(bls381_hash::probe::run(id, count));
        }
        // Field-primitive CU probes: per-op cost of mont_mul / mul2.
        22 => {
            expect_len(payload, 8)?;
            let count = u64::from_le_bytes(payload.try_into().unwrap());
            black_box(bls381_hash::probe::mont_mul_loop(count));
        }
        23 => {
            expect_len(payload, 8)?;
            let count = u64::from_le_bytes(payload.try_into().unwrap());
            black_box(bls381_hash::probe::mul2_loop(count));
        }
        _ => return Err(ProgramError::InvalidInstructionData),
    }

    Ok(())
}

type FieldG2 = <G2Projective as MapToCurve>::Field;

fn hash_to_field_g2(payload: &[u8]) -> [FieldG2; 2] {
    let mut u = [FieldG2::default(); 2];
    FieldG2::hash_to_field::<ExpandMsgXmd<Sha256>>(payload, G2_RO, &mut u);
    u
}

fn expect_len(payload: &[u8], len: usize) -> Result<(), ProgramError> {
    if payload.len() != len {
        return Err(ProgramError::InvalidInstructionData);
    }
    Ok(())
}

fn return_status(rc: u64, out: &[u8]) {
    let mut data = [0u8; 1 + GT];
    data[0] = rc as u8;
    data[1..1 + out.len()].copy_from_slice(out);
    set_return_data(&data[..1 + out.len()]);
}
