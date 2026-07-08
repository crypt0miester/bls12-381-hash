//! Quadratic extension field Fp2 = Fp[u]/(u^2 + 1), built on the Fp layer.

use solana_program_error::ProgramError;

use crate::consts_g1::R;
use crate::fp::{
    add_mod, from_mont, is_zero, mont_mul, neg_mod, sub_mod, to_mont, wit48, Fp,
};

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
