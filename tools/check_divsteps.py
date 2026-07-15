#!/usr/bin/env python3
"""Model and bound-check for fp.rs::inv_divsteps (Bernstein-Yang safegcd).

Mirrors the Rust implementation lane for lane: 13 signed 30-bit lanes,
branchless 30-step divstep batches, exact-shift f/g updates and mod-p
d/e updates. Verifies against pow(x, p-2, p) over edges and randoms,
tracks the worst accumulator magnitudes against the i64 budget, and
records the iteration bound argument.

Bound: the 2024 convex-hull analysis (Bernstein, Chen, Harrison,
Maxwell, Wang, Wuille, Yang; hull computations in sipa/safegcd-bounds)
bounds delta = 1 divsteps for 0 <= g <= f <= M, M >= 157, by
max(floor((2455 log2(M) + 1402)/1736) * 2,
    floor((2455 log2(M) + 1676)/1736) * 2 - 1):
1078 for 381 bits, so at most 36 batches run before the g == 0 exit;
the 37-batch cap keeps the older floor((49b + 57)/17) = 1101 bound's
margin. The half-delta variant (delta = 1/2, bound 878) was measured
and rejected: its smaller average w-round limit costs more than the
saved steps (see BENCHMARKS.md).
"""

import random

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
M30 = (1 << 30) - 1
PINV30 = pow(P, -1, 1 << 30)
BATCHES = 37

assert max((2455 * 381 + 1402) // 1736 * 2, (2455 * 381 + 1676) // 1736 * 2 - 1) == 1078
assert (49 * 381 + 57) // 17 == 1101
assert BATCHES * 30 >= 1101

def to_lanes(x):
    # signed-30 representation: lanes 0..11 in [0, 2^30), lane 12 signed
    lanes = [(x >> (30 * i)) & M30 for i in range(12)]
    lanes.append(x >> 360)
    return lanes

def from_lanes(l):
    return sum(v << (30 * i) for i, v in enumerate(l))

def i64(x):
    x &= (1 << 64) - 1
    return x - (1 << 64) if x >= (1 << 63) else x

class Stats:
    max_abs_cd = 0
    max_abs_md = 0
    max_row_sum = 0

def divsteps30(eta, f0, g0):
    """30 divsteps on the low lanes, batched: one multiple of f cancels up
    to min(eta+1, remaining, 6) low bits of g per round, using the odd-f
    identity f*(f^2-2) == -f^-1 mod 2^6, and the cancelled bits shift out
    in one go. Returns the transition matrix scaled by 2^30 (the f row
    doubles per consumed step instead of g halving). eta = -delta;
    reference semantics per step (checked below against the paper):
      delta > 0 and g odd: (delta, f, g) <- (1 - delta, g, (g - f)/2)
      g odd:               (delta, f, g) <- (1 + delta, f, (g + f)/2)
      g even:              (delta, f, g) <- (1 + delta, f, g/2)
    An odd add at consumed step j needs delta + j - 1 <= 0, hence the
    eta + 1 cap on a round.
    """
    u, v, q, r = 1, 0, 0, 1
    f, g = f0, g0
    i = 30
    while True:
        assert f & 1 == 1
        # strip trailing zeros of g one at a time, up to the step budget
        while g & 1 == 0:
            g >>= 1
            u <<= 1
            v <<= 1
            eta -= 1
            i -= 1
            if i == 0:
                assert max(abs(u) + abs(v), abs(q) + abs(r)) <= 1 << 30
                return eta, u, v, q, r
        if eta < 0:
            eta = -eta
            f, g = g, -f
            u, q = q, -u
            v, r = r, -v
        # one w-round: the first consumed step is the odd add, the rest
        # are the even halvings of the freshly cancelled bits
        limit = min(eta + 1, i, 6)
        m = (1 << limit) - 1
        w = (g * (f * (f * f - 2))) & m
        assert (w * f + g) & m == 0 and w & 1 == 1
        g += w * f
        q += w * u
        r += w * v
        g >>= limit
        u <<= limit
        v <<= limit
        eta -= limit
        i -= limit
        Stats.max_row_sum = max(Stats.max_row_sum, abs(u) + abs(v), abs(q) + abs(r))
        if i == 0:
            assert max(abs(u) + abs(v), abs(q) + abs(r)) <= 1 << 30
            return eta, u, v, q, r

def divstep_reference(delta, f, g):
    if delta > 0 and g & 1:
        return 1 - delta, g, (g - f) >> 1
    if g & 1:
        return 1 + delta, f, (g + f) >> 1
    return 1 + delta, f, g >> 1

def check_divsteps_vs_reference(trials=2000):
    for _ in range(trials):
        f0 = random.getrandbits(30) | 1
        g0 = random.getrandbits(30)
        delta = random.randrange(-40, 41)
        eta, u, v, q, r = divsteps30(-delta, f0, g0)
        fr, gr, dr = f0, g0, delta
        for _ in range(30):
            dr, fr, gr = divstep_reference(dr, fr, gr)
        # matrix applied to the low words reproduces the reference f, g
        assert u * f0 + v * g0 == fr << 30, "f row mismatch"
        assert q * f0 + r * g0 == gr << 30, "g row mismatch"
        assert eta == -dr, "eta mismatch"

def update_fg(f, g, fout, u, v, q, r):
    """Two row passes, mirroring the Rust register-pressure split: the f
    row lands in fout (f itself stays for the g pass), the g row updates
    in place with the shifted-store pattern."""
    cf = u * f[0] + v * g[0]
    assert cf & M30 == 0
    cf >>= 30
    for i in range(1, 13):
        cf += u * f[i] + v * g[i]
        Stats.max_abs_cd = max(Stats.max_abs_cd, abs(cf))
        fout[i - 1] = cf & M30
        cf >>= 30
    fout[12] = cf
    cg = q * f[0] + r * g[0]
    assert cg & M30 == 0
    cg >>= 30
    nonzero = 0
    for i in range(1, 13):
        cg += q * f[i] + r * g[i]
        Stats.max_abs_cd = max(Stats.max_abs_cd, abs(cg))
        g[i - 1] = cg & M30
        nonzero |= cg & M30
        cg >>= 30
    g[12] = cg
    return (nonzero | cg) == 0

def update_de(d, e, dout, u, v, q, r):
    """Two row passes like update_fg; the head quotients need the old
    d[0], e[0], so both are fixed before either pass stores."""
    sd = -1 if d[12] < 0 else 0
    se = -1 if e[12] < 0 else 0
    md = (u & sd) + (v & se)
    me = (q & sd) + (r & se)
    cd = u * d[0] + v * e[0]
    ce = q * d[0] + r * e[0]
    md -= (PINV30 * cd + md) & M30
    me -= (PINV30 * ce + me) & M30
    Stats.max_abs_md = max(Stats.max_abs_md, abs(md), abs(me))
    cd += P30[0] * md
    ce += P30[0] * me
    assert cd & M30 == 0 and ce & M30 == 0
    cd >>= 30
    ce >>= 30
    for i in range(1, 13):
        cd += u * d[i] + v * e[i] + P30[i] * md
        Stats.max_abs_cd = max(Stats.max_abs_cd, abs(cd))
        dout[i - 1] = cd & M30
        cd >>= 30
    dout[12] = cd
    for i in range(1, 13):
        ce += q * d[i] + r * e[i] + P30[i] * me
        Stats.max_abs_cd = max(Stats.max_abs_cd, abs(ce))
        e[i - 1] = ce & M30
        ce >>= 30
    e[12] = ce
    # range invariant the Rust comments claim: d, e stay in (-2p, p)
    assert -2 * P < from_lanes(dout) < P
    assert -2 * P < from_lanes(e) < P

P30 = to_lanes(P)

def inv_divsteps(a):
    assert 0 < a < P
    d, dout = to_lanes(0), to_lanes(0)
    e = to_lanes(1)
    f, fout = to_lanes(P), to_lanes(0)
    g = to_lanes(a)
    eta = -1
    done_at = None
    for batch in range(BATCHES):
        eta, u, v, q, r = divsteps30(eta, f[0], g[0])
        done = update_fg(f, g, fout, u, v, q, r)
        update_de(d, e, dout, u, v, q, r)
        f, fout = fout, f
        d, dout = dout, d
        # once g == 0 every further batch is the identity on f, d, e
        # (matrix [[2^30, 0], [0, 1]], cancelled by the shared /2^30)
        if done:
            done_at = (batch + 1) * 30
            break
    assert from_lanes(g) == 0, "g nonzero after the batch cap"
    fv = from_lanes(f)
    assert fv in (1, -1), f"gcd wandered: f = {fv}"
    dv = from_lanes(d) * fv
    assert -2 * P < dv < 2 * P
    return dv % P, done_at

def main():
    random.seed(0x9E3779B97F4A7C15)
    check_divsteps_vs_reference()
    print("divsteps30 == 30x reference divstep on 2000 random (delta, f, g)")

    worst_done = 0
    cases = [1, 2, 3, P - 1, P - 2, (P + 1) // 2, M30, 1 << 30, 1 << 360,
             (1 << 381) % P, pow(2, 390, P), pow(2, 390 * 2, P)]
    cases += [random.randrange(1, P) for _ in range(4000)]
    # adversarial-ish: long carry / sparse patterns
    cases += [(1 << k) % P for k in range(1, 381, 7)]
    cases += [((1 << k) - 1) % P for k in range(300, 381, 3)]
    for a in cases:
        if a == 0:
            continue
        inv, done = inv_divsteps(a)
        assert inv == pow(a, P - 2, P), f"wrong inverse for {a:#x}"
        worst_done = max(worst_done, done or 0)
    print(f"{len(cases)} inverses match pow(a, p-2, p)")
    print(f"worst observed divsteps to g == 0: {worst_done} (cap {BATCHES * 30}, bound 1078)")
    print(f"max |row sum| in divsteps30:       2^{Stats.max_row_sum.bit_length() - 1} <= 2^30")
    print(f"max |accumulator| in updates:      {Stats.max_abs_cd.bit_length()} bits (i64 budget 63)")
    print(f"max |md| in update_de:             {Stats.max_abs_md.bit_length()} bits (< 2^31)")

if __name__ == "__main__":
    main()
