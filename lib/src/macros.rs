//! Token expansion helpers for the generated field arithmetic

/// Multiply-accumulate a zipped run of lanes: sum += x[i] * y[j] for each
/// index pair. Pure token expansion, so the compiled code is identical to
/// writing the run out by hand; the two index lists must have the same
/// length or the expansion fails to compile.
macro_rules! dot {
    ($sum:ident, $x:ident $($i:literal)*, $y:ident $($j:literal)*) => {
        $( $sum += $x[$i] * $y[$j]; )*
    };
}

/// Pick the quotient lane m[k] that zeroes the low 30 bits, then shift
macro_rules! quotient {
    ($sum:ident, $m:ident $k:literal) => {
        $m[$k] = $sum.wrapping_mul(INV30) & MASK30;
        $sum = ($sum + $m[$k] * P30[0]) >> 30;
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

pub(crate) use dot;
pub(crate) use lane;
pub(crate) use quotient;
pub(crate) use stage;
