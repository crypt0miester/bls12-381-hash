#!/usr/bin/env python3
"""Experiment record: 60-divstep composed updates. VERDICT: rejected.

Measured 2026-07-14 on SBF v3 (probe tag 24 id 26, since removed):
34,590.6 CU per inverse against 34,740.0 for the shipping 30-step form,
a 149 CU wash. Composition conserves the product count, so the win had
to come from halved pass overhead, and the fused row pass (four split
entries, two quotient limbs, four sliding operands, one accumulator)
exceeds the ten sBPF registers; the spills give the overhead back. The
algorithm below is correct (verified against Fermat) and kept as the
record.

Idea from the safegcd 2024 follow-up talk: run two 30-step w-round batches
back to back (refreshing the low words in between with a two-lane partial
update), compose the transition matrices into 60-bit entries, and apply
them to f, g, d, e in ONE pass per row with balanced 30-bit entry splits
and a two-limb quotient. Halves the number of full-width passes; product
count stays the same, the hoped saving is per-pass overhead.

Columns must fit i64. Balanced splits give, per column:
  |lo terms|  <= 2*2^29*2^30 = 2^60
  |hi terms|  <= 2*(2^30+1)*2^30 ~ 2^61
  |quotients| <= 2*2^29*2^30 = 2^60   (balanced md limbs)
  total ~ 2^61.9 + carry, inside the 2^63 budget; this script tracks the
  actual worst case. The d/e range needs one conditional +-p at the end
  (the sign trick does not transfer to signed quotients); the script
  tracks the worst |d|, |e| against the 2p target.
"""

import random

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
M30 = (1 << 30) - 1
PINV30 = pow(P, -1, 1 << 30)
BATCHES60 = 19  # 19*60 = 1140 >= 1101 (and >= 1078)

def to_lanes(x):
    lanes = [(x >> (30 * i)) & M30 for i in range(12)]
    lanes.append(x >> 360)
    return lanes

def from_lanes(l):
    return sum(v << (30 * i) for i, v in enumerate(l))

P30 = to_lanes(P)

class Stats:
    max_abs_col = 0
    max_abs_ratio_p = 0.0
    max_entry = 0

def divsteps30(eta, f0, g0):
    u, v, q, r = 1, 0, 0, 1
    f, g = f0, g0
    i = 30
    while True:
        while g & 1 == 0:
            g >>= 1
            u <<= 1
            v <<= 1
            eta -= 1
            i -= 1
            if i == 0:
                return eta, u, v, q, r
        if eta < 0:
            eta = -eta
            f, g = g, -f
            u, q = q, -u
            v, r = r, -v
        limit = min(eta + 1, i, 6)
        m = (1 << limit) - 1
        w = (g * (f * (f * f - 2))) & m
        g += w * f
        q += w * u
        r += w * v
        g >>= limit
        u <<= limit
        v <<= limit
        eta -= limit
        i -= limit
        if i == 0:
            return eta, u, v, q, r

def divsteps60(eta, f, g):
    """Two 30-step batches with a two-lane low-word refresh in between,
    composed into one 60-bit matrix."""
    eta, u1, v1, q1, r1 = divsteps30(eta, f[0], g[0])
    # low 30 bits of the stepped f, g from lanes 0..1 (exact: the matrix
    # zeroes the low 30 bits, lane products stay under 2^61)
    sf0 = u1 * f[0] + v1 * g[0]
    sg0 = q1 * f[0] + r1 * g[0]
    assert sf0 & M30 == 0 and sg0 & M30 == 0
    f2 = ((sf0 >> 30) + u1 * f[1] + v1 * g[1]) & M30
    g2 = ((sg0 >> 30) + q1 * f[1] + r1 * g[1]) & M30
    eta, u2, v2, q2, r2 = divsteps30(eta, f2, g2)
    u = u2 * u1 + v2 * q1
    v = u2 * v1 + v2 * r1
    q = q2 * u1 + r2 * q1
    r = q2 * v1 + r2 * r1
    Stats.max_entry = max(Stats.max_entry, abs(u) + abs(v), abs(q) + abs(r))
    assert abs(u) + abs(v) <= 1 << 60 and abs(q) + abs(r) <= 1 << 60
    return eta, u, v, q, r

def bal_split(x):
    """x = hi*2^30 + lo with lo in [-2^29, 2^29)."""
    lo = ((x + (1 << 29)) & M30) - (1 << 29)
    hi = (x - lo) >> 30
    assert hi * (1 << 30) + lo == x
    return hi, lo

def row60(a, b, x, y, quot=False):
    """One output row t = (x*a + y*b [+ mp*P]) / 2^60 over 13 lanes with
    balanced splits; returns lanes. Tracks worst column magnitude."""
    x_hi, x_lo = bal_split(x)
    y_hi, y_lo = bal_split(y)
    # head: columns 0 and 1 must vanish mod 2^60 (after quotient limbs)
    col0 = x_lo * a[0] + y_lo * b[0]
    m0 = 0
    if quot:
        m0 = ((-col0 * PINV30) + (1 << 29)) % (1 << 30) - (1 << 29)
        col0 += m0 * P30[0]
    assert col0 & M30 == 0
    c = col0 >> 30
    col1 = c + x_lo * a[1] + y_lo * b[1] + x_hi * a[0] + y_hi * b[0]
    m1 = 0
    if quot:
        col1 += m0 * P30[1]
        m1 = ((-col1 * PINV30) + (1 << 29)) % (1 << 30) - (1 << 29)
        col1 += m1 * P30[0]
    assert col1 & M30 == 0
    c = col1 >> 30
    out = [0] * 13
    for i in range(2, 13):
        col = c + x_lo * a[i] + y_lo * b[i] + x_hi * a[i - 1] + y_hi * b[i - 1]
        if quot:
            col += m0 * P30[i] + m1 * P30[i - 1]
        Stats.max_abs_col = max(Stats.max_abs_col, abs(col))
        out[i - 2] = col & M30
        c = col >> 30
    # tail columns 13 and 14 (only hi terms and quotient tails)
    col = c + x_hi * a[12] + y_hi * b[12]
    if quot:
        col += m0 * 0 + m1 * P30[12]
    Stats.max_abs_col = max(Stats.max_abs_col, abs(col))
    out[11] = col & M30
    out[12] = col >> 30
    return out

def recenter(t, sign):
    """t += sign*p as a lane pass with signed carries, like the Rust."""
    c = 0
    for i in range(12):
        lane = t[i] + sign * P30[i] + c
        t[i] = lane & M30
        c = lane >> 30
    t[12] = t[12] + sign * P30[12] + c

def inv_divsteps60(a):
    assert 0 < a < P
    d = to_lanes(0)
    e = to_lanes(1)
    f = to_lanes(P)
    g = to_lanes(a)
    eta = -1
    done = False
    for _ in range(BATCHES60):
        eta, u, v, q, r = divsteps60(eta, f, g)
        nf = row60(f, g, u, v)
        ng = row60(f, g, q, r)
        nd = row60(d, e, u, v, quot=True)
        ne = row60(d, e, q, r, quot=True)
        # conditional +-p recentering replaces the sign trick: trigger on
        # the top lane alone (the escape window above 2p is covered by the
        # divergence analysis, worst |d| stays under ~2.6p)
        for t in (nd, ne):
            if t[12] > 2 * P30[12]:
                recenter(t, -1)
            elif t[12] < -2 * P30[12]:
                recenter(t, 1)
        f, g, d, e = nf, ng, nd, ne
        Stats.max_abs_ratio_p = max(
            Stats.max_abs_ratio_p, abs(from_lanes(d)) / P, abs(from_lanes(e)) / P
        )
        assert Stats.max_abs_ratio_p < 3.0
        if from_lanes(g) == 0:
            done = True
            break
    assert done, "g nonzero after cap"
    fv = from_lanes(f)
    assert fv in (1, -1)
    return (from_lanes(d) * fv) % P

def main():
    random.seed(0x9E3779B97F4A7C15)
    cases = [1, 2, 3, P - 1, P - 2, (P + 1) // 2, M30, 1 << 30, 1 << 360]
    cases += [random.randrange(1, P) for _ in range(2000)]
    for a in cases:
        assert inv_divsteps60(a) == pow(a, P - 2, P), f"wrong inverse {a:#x}"
    print(f"{len(cases)} inverses match pow(a, p-2, p)")
    print(f"max row-sum of 60-bit entries: 2^{Stats.max_entry.bit_length() - 1}")
    print(f"max |column| in row60:         {Stats.max_abs_col.bit_length()} bits (budget 63)")
    print(f"max |d|,|e| in units of p:     {Stats.max_abs_ratio_p:.3f} (lanes hold ~4)")

if __name__ == "__main__":
    main()
