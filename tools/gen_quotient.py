# Straight-line body for quotient_matches: lanes of x (8 limbs, 18 lanes),
# lanes of q (3 limbs, 5 lanes), then 18 columns of u + q * P30.
def lane_expr(src, j, nlimbs):
    bit = 30 * j; l, s = divmod(bit, 64)
    if l >= nlimbs:
        return "0"
    e = f"({src}[{l}] >> {s})" if s else f"{src}[{l}]"
    if s > 34 and l + 1 < nlimbs:
        e = f"({e} | ({src}[{l+1}] << {64 - s}))"
    return f"{e} & MASK30"
out = []
out.append("    let x = [")
for j in range(18):
    out.append(f"        {lane_expr('limbs', j, 8)},")
out.append("    ];")
out.append("    let q = [")
for j in range(5):
    out.append(f"        {lane_expr('ql', j, 3)},")
out.append("    ];")
out.append("    let u = split30(u);")
out.append("    let mut col = 0u64;")
out.append("    let mut diff = 0u64;")
for j in range(18):
    ks = [k for k in range(5) if 0 <= j - k < 13]
    if j < 13:
        out.append(f"    col = col.wrapping_add(u[{j}]);")
    qi = " ".join(str(k) for k in ks); pi = " ".join(str(j - k) for k in ks)
    if ks: out.append(f"    dot!(col, q {qi}, P30 {pi});")
    out.append(f"    diff |= (col & MASK30) ^ x[{j}];")
    out.append(f"    col >>= 30;")
out.append("    (diff | col) == 0")
print("\n".join(out))
