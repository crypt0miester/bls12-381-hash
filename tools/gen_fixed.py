# Straight-line body for a fixed-multiplier Montgomery product:
# sum_i a_i T[i] over 13 lanes, T[i] = K 2^(30 i - 330) mod p, then two
# Montgomery quotient steps divide by 2^60. Emits a macro_rules! body.
lines = []
lines.append("macro_rules! fixed_mul {")
lines.append("    ($t:expr, $a:ident) => {{")
lines.append("        let a = split30($a);")
lines.append("        let mut r = [0u64; 13];")
lines.append("        let mut col = 0u64;")
def dot(j):
    return [f"        col = col.wrapping_add(a[{i}].wrapping_mul($t[{i}][{j}]));" for i in range(13)]
# col 0
lines += dot(0)
lines.append("        let m0 = col.wrapping_mul(INV30) & MASK30;")
lines.append("        col = col.wrapping_add(m0.wrapping_mul(P30[0])) >> 30;")
# col 1
lines += dot(1)
lines.append("        col = col.wrapping_add(m0.wrapping_mul(P30[1]));")
lines.append("        let m1 = col.wrapping_mul(INV30) & MASK30;")
lines.append("        col = col.wrapping_add(m1.wrapping_mul(P30[0])) >> 30;")
for j in range(2, 13):
    lines += dot(j)
    lines.append(f"        col = col.wrapping_add(m0.wrapping_mul(P30[{j}])).wrapping_add(m1.wrapping_mul(P30[{j-1}]));")
    lines.append(f"        r[{j-2}] = col & MASK30;")
    lines.append("        col >>= 30;")
lines.append("        col = col.wrapping_add(m1.wrapping_mul(P30[12]));")
lines.append("        r[11] = col & MASK30;")
lines.append("        r[12] = col >> 30;")
lines.append("        pack30(&r)")
lines.append("    }};")
lines.append("}")
print("\n".join(lines))
