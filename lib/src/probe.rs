//! Per function CU probes for the bench harness; not part of the API.
//! Every loop carries a data dependency so the optimizer cannot hoist or
//! collapse the body, and each id isolates exactly one function.

// Ids 6 and 22 fold at compile time (negation and the byte round trip are
// involutions), so their numbers are floors, not costs.
use crate::consts_g1::{R, R2};
use crate::fp::{
    add_mod, be_to_limbs, from_mont, inv_divsteps, limbs_to_be, mont_mul, mont_sqr, neg_mod,
    sub_mod, to_mont, wit48, Fp,
};
use crate::fp2::{add2, from_mont2, mul2, sq2, sub2, to_mont2, wit96, Fp2};
use crate::g1::{expand_message_xmd, gx_at, iso11_adapted, mul_by_xi};
use crate::g2::{expand_message_xmd_g2, gx2_at, iso3_adapted, mul_by_a2i, mul_by_xi2};

pub fn mont_mul_loop(count: u64) -> u64 {
    let mut a = R;
    for _ in 0..count {
        a = mont_mul(&a, &R2);
    }
    a[0]
}

pub fn mul2_loop(count: u64) -> u64 {
    let mut a = Fp2 { c0: R, c1: R2 };
    let b = Fp2 { c0: R2, c1: R };
    for _ in 0..count {
        a = mul2(&a, &b);
    }
    a.c0[0]
}

fn fp_loop(count: u64, op: impl Fn(&Fp) -> Fp) -> u64 {
    let mut a = R;
    for _ in 0..count {
        a = op(&a);
    }
    a[0]
}

fn fp2_loop(count: u64, op: impl Fn(&Fp2) -> Fp2) -> u64 {
    let mut a = Fp2 { c0: R, c1: R2 };
    for _ in 0..count {
        a = op(&a);
    }
    a.c0[0]
}

pub fn run(id: u8, count: u64) -> u64 {
    let b2 = Fp2 { c0: R2, c1: R };
    match id {
        0 => mont_mul_loop(count),
        1 => fp_loop(count, mont_sqr),
        2 => fp_loop(count, from_mont),
        3 => fp_loop(count, to_mont),
        4 => fp_loop(count, |a| add_mod(a, &R2)),
        5 => fp_loop(count, |a| sub_mod(a, &R2)),
        6 => fp_loop(count, neg_mod),
        7 => mul2_loop(count),
        8 => fp2_loop(count, sq2),
        9 => fp2_loop(count, |a| add2(a, &b2)),
        10 => fp2_loop(count, |a| sub2(a, &b2)),
        11 => fp2_loop(count, mul_by_xi2),
        12 => fp2_loop(count, mul_by_a2i),
        13 => fp2_loop(count, from_mont2),
        14 => fp2_loop(count, to_mont2),
        15 => fp_loop(count, |x| {
            let (a, b, c, d) = iso11_adapted(x);
            add_mod(&add_mod(&a, &b), &add_mod(&c, &d))
        }),
        16 => fp2_loop(count, |x| {
            let (a, b, c, d) = iso3_adapted(x);
            add2(&add2(&a, &b), &add2(&c, &d))
        }),
        17 => fp_loop(count, gx_at),
        18 => fp2_loop(count, gx2_at),
        19 => fp_loop(count, mul_by_xi),
        20 => {
            // wit48: parse and range check a 48 byte canonical witness
            let mut acc = [0u64; 6];
            let mut bytes = [0u8; 48];
            for i in 0..count {
                bytes[47] = acc[0] as u8 ^ i as u8;
                acc = add_mod(&acc, &wit48(&bytes).unwrap());
            }
            acc[0]
        }
        21 => {
            let mut acc = [0u64; 6];
            let mut bytes = [0u8; 96];
            for i in 0..count {
                bytes[47] = acc[0] as u8 ^ i as u8;
                bytes[95] = bytes[47].wrapping_add(1);
                let v = wit96(&bytes).unwrap();
                acc = add_mod(&acc, &v.c0);
            }
            acc[0]
        }
        22 => {
            // serialization round trip
            let mut a = R;
            for _ in 0..count {
                a = be_to_limbs(&limbs_to_be(&a));
                a[0] ^= 1;
            }
            a[0]
        }
        23 => {
            let mut acc = 7u64;
            for _ in 0..count {
                let b = expand_message_xmd(crate::dst::G1_RO, &acc.to_le_bytes());
                acc = acc.wrapping_add(b[0][0] as u64);
            }
            acc
        }
        24 => {
            let mut acc = 7u64;
            for _ in 0..count {
                let b = expand_message_xmd_g2(crate::dst::G2_RO, &acc.to_le_bytes());
                acc = acc.wrapping_add(b[0][0] as u64);
            }
            acc
        }
        25 => {
            // divsteps inverse; alternates a and a^-1, loop-carried
            let mut a = R;
            for _ in 0..count {
                a = inv_divsteps(&a).unwrap();
            }
            a[0]
        }
        _ => 0,
    }
}
