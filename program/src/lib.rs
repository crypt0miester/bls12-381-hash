#![allow(unexpected_cfgs)]

mod g1_consts;
pub mod g1_msig;
mod g2_consts;
pub mod g2_msig;

use core::hint::black_box;

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

const DST_G2: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_";
const DST_G1: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_POP_";

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
                payload, DST_G2,
            );
            set_return_data(&G2Affine::from(p).to_compressed());
        }
        1 => {
            let p = <G1Projective as HashToCurve<ExpandMsgXmd<Sha256>>>::hash_to_curve(
                payload, DST_G1,
            );
            set_return_data(&G1Affine::from(p).to_compressed());
        }
        2 => {
            let u = hash_to_field_g2(payload);
            black_box(&u);
        }
        3 => {
            let u = hash_to_field_g2(payload);
            let p = G2Projective::map_to_curve(&u[0]) + &G2Projective::map_to_curve(&u[1]);
            black_box(&p);
        }
        4 => {
            let u = hash_to_field_g2(payload);
            let p = (G2Projective::map_to_curve(&u[0]) + &G2Projective::map_to_curve(&u[1]))
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
            let out = g1_msig::run(tag - 30, payload)?;
            set_return_data(&out);
        }
        34 => {
            expect_len(payload, 8)?;
            let count = u64::from_le_bytes(payload.try_into().unwrap());
            let acc = g1_msig::mont_mul_bench(count);
            set_return_data(&acc.to_le_bytes());
        }
        // End-to-end min-pk vote verify: witness hash_to_G2, subtract the
        // absentees from the stored committee aggregate, check the pairing
        // equation e(effective_pk, H) == e(g1_gen, sig).
        51 => {
            let absent = payload[0] as usize;
            let agg_end = 1 + G1_POINT;
            let gen_end = agg_end + G1_POINT;
            let sig_end = gen_end + G2_COMPRESSED;
            let abs_end = sig_end + 48 * absent;
            if payload.len() < abs_end {
                return Err(ProgramError::InvalidInstructionData);
            }

            let mut effective: [u8; G1_POINT] = payload[1..agg_end].try_into().unwrap();
            for i in 0..absent {
                let compressed = &payload[sig_end + i * 48..sig_end + (i + 1) * 48];
                let mut member = [0u8; G1_POINT];
                let rc = unsafe {
                    sys::sol_curve_decompress(
                        BLS12_381_G1_BE,
                        compressed.as_ptr(),
                        member.as_mut_ptr(),
                    )
                };
                if rc != 0 {
                    return Err(ProgramError::InvalidInstructionData);
                }
                let mut out = [0u8; G1_POINT];
                let rc = unsafe {
                    sys::sol_curve_group_op(
                        BLS12_381_G1_BE,
                        OP_SUB,
                        effective.as_ptr(),
                        member.as_ptr(),
                        out.as_mut_ptr(),
                    )
                };
                if rc != 0 {
                    return Err(ProgramError::InvalidInstructionData);
                }
                effective = out;
            }

            let mut sig = [0u8; G2_POINT];
            let rc = unsafe {
                sys::sol_curve_decompress(
                    BLS12_381_G2_BE,
                    payload[gen_end..sig_end].as_ptr(),
                    sig.as_mut_ptr(),
                )
            };
            if rc != 0 {
                return Err(ProgramError::InvalidInstructionData);
            }

            let hash = g2_msig::run_witnessed(&payload[abs_end..])?;

            let mut gt_left = [0u8; GT];
            let rc = unsafe {
                sys::sol_curve_pairing_map(
                    BLS12_381_PAIRING_BE,
                    1,
                    effective.as_ptr(),
                    hash.as_ptr(),
                    gt_left.as_mut_ptr(),
                )
            };
            if rc != 0 {
                return Err(ProgramError::InvalidInstructionData);
            }
            let mut gt_right = [0u8; GT];
            let rc = unsafe {
                sys::sol_curve_pairing_map(
                    BLS12_381_PAIRING_BE,
                    1,
                    payload[agg_end..gen_end].as_ptr(),
                    sig.as_ptr(),
                    gt_right.as_mut_ptr(),
                )
            };
            if rc != 0 || gt_left != gt_right {
                return Err(ProgramError::InvalidInstructionData);
            }
            set_return_data(&[1]);
        }
        40 => {
            let out = g1_msig::run_witnessed(payload)?;
            set_return_data(&out);
        }
        41 => {
            let out = g2_msig::run_witnessed(payload)?;
            set_return_data(&out);
        }
        35 => {
            expect_len(payload, 9)?;
            let count = u64::from_le_bytes(payload[1..].try_into().unwrap());
            let acc = g1_msig::mul_bench(payload[0], count);
            set_return_data(&acc.to_le_bytes());
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
        _ => return Err(ProgramError::InvalidInstructionData),
    }

    Ok(())
}

type FieldG2 = <G2Projective as MapToCurve>::Field;

fn hash_to_field_g2(payload: &[u8]) -> [FieldG2; 2] {
    let mut u = [FieldG2::default(); 2];
    FieldG2::hash_to_field::<ExpandMsgXmd<Sha256>>(payload, DST_G2, &mut u);
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
