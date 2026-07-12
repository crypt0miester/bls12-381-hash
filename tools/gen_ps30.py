#!/usr/bin/env python3
# Emits the straight line bodies of mont_mul, mont_sqr and from_mont in
# lib/src/fp.rs. The bodies are spelled out because the sbpf backend leaves
# loops rolled.
#
# ps30: thirteen 30 bit lanes, R = 2^390. Against ps29 that is one lane
# fewer, cutting the multiply from 392 to 338 products, but the tall middle
# columns can now sum past 2^64. The generator tracks the exact worst case
# of the running column sum (modulus lanes exact, operand lanes bounded by
# the 2^384 operand contract) and, where the full column would overflow,
# banks the high bits into a spill register between the operand run and the
# modulus run. Both halves are asserted to fit; so is every carry.
# Run and paste over the three functions:
#   python3 tools/gen_ps30.py
P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
LANES = 13
WIDTH = 30
MASK = (1 << WIDTH) - 1
U64 = (1 << 64) - 1
OPERAND = (1 << 384) - 1  # mont_mul digests any value below 2^384

inverse = pow(P, -1, 1 << WIDTH)
NEG_INV = (1 << WIDTH) - inverse
assert (P * NEG_INV + 1) % (1 << WIDTH) == 0
assert OPERAND * OPERAND < P << (LANES * WIDTH)  # a*b < p*R keeps the output below 2p


def split(x):
    return [(x >> (WIDTH * i)) & MASK for i in range(LANES)]


P_LANE = split(P)
A_LANE = split(OPERAND)  # per lane operand ceiling; top lane is 2^24 - 1
M_LANE = [MASK] * LANES  # quotient lanes are masked to WIDTH bits
assert all(l < (1 << 31) for l in P_LANE)  # valid positive imm32 operands


def dot_line(x, y, k, first, last):
    if first > last:
        return []
    left = " ".join(str(i) for i in range(first, last + 1))
    right = " ".join(str(k - i) for i in range(first, last + 1))
    return [f"    dot!(sum, {x} {left}, {y} {right});"]


SPILL = ["    let spill = sum >> 30;", "    sum &= MASK30;"]


class Column:
    def __init__(self, k, lines, ab_max, mp_max):
        self.k = k
        self.lines = lines  # operand-run source lines
        self.ab_max = ab_max  # worst case sum of the operand run
        self.mp_max = mp_max  # worst case sum of the modulus run (quotient included)


def emit(columns):
    out = []
    spills = set()
    carry = 0
    for c in columns:
        out.append(f"    // column {c.k}")
        out += c.lines
        low = c.k < LANES
        mp_lines = (
            dot_line("m", "P30", c.k, 0, c.k - 1)
            if low
            else dot_line("m", "P30", c.k, c.k - LANES + 1, LANES - 1)
        )
        ab_total = carry + c.ab_max
        assert ab_total <= U64, f"column {c.k}: operand run alone overflows"
        if ab_total + c.mp_max > U64:
            out += SPILL
            spills.add(c.k)
            spilled = ab_total >> WIDTH
            total = MASK + c.mp_max
            assert total <= U64, f"column {c.k}: modulus run overflows after spill"
        else:
            spilled = 0
            total = ab_total + c.mp_max
        out += mp_lines
        if low:
            out.append(f"    quotient!(sum, m {c.k});")
        else:
            out.append(f"    lane!(sum, r {c.k - LANES});")
        if spilled:
            out.append("    sum += spill;")
        carry = (total >> WIDTH) + spilled
    assert carry <= MASK, "final lane exceeds one lane"
    return out, spills


def simulate(spills, terms_for, value):
    # Execute exactly what the emitted code does, u64 overflow asserted
    m = [0] * LANES
    r = [0] * LANES
    s = 0

    def add(v):
        nonlocal s
        s += v
        assert s <= U64

    for k in range(2 * LANES - 1):
        for v in terms_for(k):
            add(v)
        spill = 0
        if k in spills:
            spill, s = s >> WIDTH, s & MASK
        if k < LANES:
            for i in range(k):
                add(m[i] * P_LANE[k - i])
            m[k] = (s * NEG_INV) & MASK
            add(m[k] * P_LANE[0])
            assert s & MASK == 0
            s >>= WIDTH
        else:
            for i in range(k - LANES + 1, LANES):
                add(m[i] * P_LANE[k - i])
            r[k - LANES] = s & MASK
            s >>= WIDTH
        s += spill
        assert s <= U64
    assert s <= MASK
    r[LANES - 1] = s
    packed = sum(l << (WIDTH * i) for i, l in enumerate(r))
    out = packed - P if packed >= P else packed
    assert 0 <= out < P
    assert (out << (LANES * WIDTH)) % P == value % P
    return out


def selftest(mul_spills, sqr_spills, redc_spills):
    import random

    random.seed(29 * 30)
    edges = [0, 1, P - 1, P, P + 1, 2 * P, 5 * P, OPERAND, OPERAND - P, 1 << 383]
    samples = edges + [random.randrange(OPERAND + 1) for _ in range(300)]
    for a in samples:
        al = split(a)
        for b in random.sample(samples, 40):
            bl = split(b)
            simulate(
                mul_spills,
                lambda k: [
                    al[i] * bl[k - i]
                    for i in range(max(0, k - LANES + 1), min(k, LANES - 1) + 1)
                ],
                a * b,
            )
        simulate(
            sqr_spills,
            lambda k: [
                2 * al[i] * al[k - i]
                for i in range(max(0, k - LANES + 1), (k + 1) // 2)
            ]
            + ([al[k // 2] ** 2] if k % 2 == 0 else []),
            a * a,
        )
        simulate(redc_spills, lambda k: [al[k]] if k < LANES else [], a)


def mul_column(k):
    first, last = max(0, k - LANES + 1), min(k, LANES - 1)
    ab_max = sum(A_LANE[i] * A_LANE[k - i] for i in range(first, last + 1))
    return Column(k, dot_line("a", "b", k, first, last), ab_max, mp_max(k))


def sqr_column(k):
    first, last = max(0, k - LANES + 1), (k + 1) // 2 - 1
    ab_max = sum(2 * A_LANE[i] * A_LANE[k - i] for i in range(first, last + 1))
    lines = dot_line("twice", "a", k, first, last)
    if k % 2 == 0:
        lines.append(f"    sum += a[{k // 2}] * a[{k // 2}];")
        ab_max += A_LANE[k // 2] ** 2
    return Column(k, lines, ab_max, mp_max(k))


def redc_column(k):
    if k < LANES:
        return Column(k, [f"    sum += t[{k}];"], A_LANE[k], mp_max(k))
    return Column(k, [], 0, mp_max(k))


def mp_max(k):
    if k < LANES:
        # m_0..m_{k-1} against the mirrored modulus lanes, then the
        # quotient product m_k * p_0
        return sum(M_LANE[i] * P_LANE[k - i] for i in range(k)) + MASK * P_LANE[0]
    return sum(M_LANE[i] * P_LANE[k - i] for i in range(k - LANES + 1, LANES))


def body(signature, prelude, column_for):
    columns = [column_for(k) for k in range(2 * LANES - 1)]
    lines, spills = emit(columns)
    out = [signature, *prelude, "", "    let mut sum = 0u64;"]
    out += lines
    out.append("    debug_assert!(sum <= MASK30);")
    out.append(f"    r[{LANES - 1}] = sum;")
    out.append("    pack30(&r)")
    out.append("}")
    return "\n".join(out), spills


mul, mul_spills = body(
    "pub(crate) fn mont_mul(a: &Fp, b: &Fp) -> Fp {",
    [
        "    let a = split30(a);",
        "    let b = split30(b);",
        "    let mut m = [0u64; 13];",
        "    let mut r = [0u64; 13];",
    ],
    mul_column,
)
sqr, sqr_spills = body(
    "pub(crate) fn mont_sqr(a: &Fp) -> Fp {",
    [
        "    let a = split30(a);",
        "    let twice = double_lanes(&a);",
        "    let mut m = [0u64; 13];",
        "    let mut r = [0u64; 13];",
    ],
    sqr_column,
)
redc, redc_spills = body(
    "/// Out of Montgomery form: a times R^-1 mod p, the bare reduction\npub(crate) fn from_mont(x: &Fp) -> Fp {",
    [
        "    let t = split30(x);",
        "    let mut m = [0u64; 13];",
        "    let mut r = [0u64; 13];",
    ],
    redc_column,
)
selftest(mul_spills, sqr_spills, redc_spills)
print(mul)
print()
print(sqr)
print()
print(redc)
