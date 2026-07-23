//! Token expansion helpers for the generated field arithmetic

/// Multiply-accumulate a zipped run of lanes: sum += x[i] * y[j] for each
/// index pair. Pure token expansion, so the compiled code is identical to
/// writing the run out by hand; the two index lists must have the same
/// length or the expansion fails to compile.
macro_rules! dot {
    ($sum:ident, $x:ident $($i:literal)*, $y:ident $($j:literal)*) => {
        $( $sum = $sum.wrapping_add($x[$i].wrapping_mul($y[$j])); )*
    };
}

/// Pick the quotient lane m[k] that zeroes the low 30 bits, then shift
macro_rules! quotient {
    ($sum:ident, $m:ident $k:literal) => {
        $m[$k] = $sum.wrapping_mul(INV30) & MASK30;
        $sum = $sum.wrapping_add($m[$k].wrapping_mul(P30[0])) >> 30;
    };
}

/// Emit one result lane and carry the rest
macro_rules! lane {
    ($sum:ident, $r:ident $j:literal) => {
        $r[$j] = $sum & MASK30;
        $sum >>= 30;
    };
}

/// One Knuth-adapted stage: t = w + k1 (+ x per extra term), then
/// acc = acc * t + k2. The sums skip the modular reduction: t stays below
/// 5p and acc below 2p, both valid mont_mul operands, and the multiply
/// output comes back canonical.
macro_rules! stage {
    ($acc:ident, $w:ident + $k1:expr $(, $x:ident)*; $k2:expr) => {{
        let t = add_unreduced(&$w, &$k1);
        $( let t = add_unreduced(&t, $x); )*
        $acc = add_unreduced(&mont_mul(&$acc, &t), &$k2);
    }};
}

/// Straight-line divsteps row pass: dst = (x a + y b) / 2^30 with the
/// shifted store. One matrix row per pass keeps the live set inside the
/// ten sBPF registers; expansion only, like dot!, because rolled loops
/// pay index arithmetic and loop control per limb.
macro_rules! row_limbs {
    ($dst:ident $a:ident $b:ident $c:ident, $x:ident $y:ident; $($i:literal)*) => {
        $(
            $c = $c.wrapping_add($x.wrapping_mul($a[$i])).wrapping_add($y.wrapping_mul($b[$i]));
            $dst[$i - 1] = $c & M30S;
            $c >>= 30;
        )*
    };
}

/// The g row pass: row_limbs plus the or-fold for the zero test
macro_rules! row_limbs_fold {
    ($dst:ident $a:ident $b:ident $c:ident $nonzero:ident, $x:ident $y:ident; $($i:literal)*) => {
        $(
            $c = $c.wrapping_add($x.wrapping_mul($a[$i])).wrapping_add($y.wrapping_mul($b[$i]));
            $dst[$i - 1] = $c & M30S;
            $nonzero |= $c & M30S;
            $c >>= 30;
        )*
    };
}

/// The d/e row pass: row_limbs plus the modulus quotient term realizing
/// the mod-p division by 2^30
macro_rules! row_limbs_modp {
    ($dst:ident $a:ident $b:ident $c:ident, $x:ident $y:ident $m:ident; $($i:literal)*) => {
        $(
            $c = $c.wrapping_add($x.wrapping_mul($a[$i])).wrapping_add($y.wrapping_mul($b[$i])).wrapping_add(P30S[$i].wrapping_mul($m));
            $dst[$i - 1] = $c & M30S;
            $c >>= 30;
        )*
    };
}

pub(crate) use dot;
pub(crate) use lane;
pub(crate) use quotient;
pub(crate) use row_limbs;
pub(crate) use row_limbs_fold;
pub(crate) use row_limbs_modp;
pub(crate) use stage;
