#!/usr/bin/env python3
"""Derive the Velu form of the RFC 9380 G2 3-isogeny and check it.

The iso-3 denominators are powers of t = x - x_k for the kernel point
x_k = -6 + 6i (x_den = t^2, y_den = t^3). Taylor-shifting the numerators to
t leaves X = a3 (x + 48i / t + 16 (1 + i) / t^2) and
Y = c a3 y (1 - 48i / t^2 - 32 (1 + i) / t^3) with a3 and c real. With
sigma = 4 / t every coefficient is an integer below 16:
  X = a3 (x + 12i sigma + (1 + i) sigma^2)
  Y = (c a3 / 2) y (2 - 6i sigma^2 - (1 + i) sigma^3)
This script rebuilds that from the RFC tables, checks it against the Horner
map on random points, checks the curve identity the fat layout leans on
(Y^2 - X^3 - B = F^2 (y^2 - g'(x)) for off-curve y too), checks the two x
pins, and prints the constants consts_g2.rs carries.
"""
import random
from math import comb

p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab


class F2:
    def __init__(s, a, b=0):
        s.a, s.b = a % p, b % p

    def __add__(s, o):
        o = lift(o)
        return F2(s.a + o.a, s.b + o.b)

    __radd__ = __add__

    def __sub__(s, o):
        o = lift(o)
        return F2(s.a - o.a, s.b - o.b)

    def __neg__(s):
        return F2(-s.a, -s.b)

    def __mul__(s, o):
        o = lift(o)
        return F2(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a)

    __rmul__ = __mul__

    def __eq__(s, o):
        o = lift(o)
        return s.a == o.a and s.b == o.b

    def inv(s):
        n = pow((s.a * s.a + s.b * s.b) % p, p - 2, p)
        return F2(s.a * n, -s.b * n)

    def __truediv__(s, o):
        return s * lift(o).inv()

    def __pow__(s, e):
        r, b = F2(1), s
        while e:
            if e & 1:
                r = r * b
            b, e = b * b, e >> 1
        return r


def lift(o):
    return o if isinstance(o, F2) else F2(o)


H = lambda h: int(h, 16)
# RFC 9380 appendix E.3, low to high degree
XNUM = [F2(H('5c759507e8e333ebb5b7a9a47d7ed8532c52d39fd3a042a88b58423c50ae15d5c2638e343d9c71c6238aaaaaaaa97d6'), H('5c759507e8e333ebb5b7a9a47d7ed8532c52d39fd3a042a88b58423c50ae15d5c2638e343d9c71c6238aaaaaaaa97d6')),
        F2(0, H('11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71a')),
        F2(H('11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71e'), H('8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffffe38d')),
        F2(H('171d6541fa38ccfaed6dea691f5fb614cb14b4e7f4e810aa22d6108f142b85757098e38d0f671c7188e2aaaaaaaa5ed1'))]
XDEN = [F2(0, H('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaa63')),
        F2(0xc, H('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaa9f')), F2(1)]
YNUM = [F2(H('1530477c7ab4113b59a4c18b076d11930f7da5d4a07f649bf54439d87d27e500fc8c25ebf8c92f6812cfc71c71c6d706'), H('1530477c7ab4113b59a4c18b076d11930f7da5d4a07f649bf54439d87d27e500fc8c25ebf8c92f6812cfc71c71c6d706')),
        F2(0, H('5c759507e8e333ebb5b7a9a47d7ed8532c52d39fd3a042a88b58423c50ae15d5c2638e343d9c71c6238aaaaaaaa97be')),
        F2(H('11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71c'), H('8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffffe38f')),
        F2(H('124c9ad43b6cf79bfbf7043de3811ad0761b0f37a1e26286b0e977c69aa274524e79097a56dc4bd9e1b371c71c718b10'))]
YDEN = [F2(H('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffa8fb'), H('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffa8fb')),
        F2(0, H('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffa9d3')),
        F2(0x12, H('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaa99')), F2(1)]
A, B, BE = F2(0, 240), F2(1012, 1012), F2(4, 4)
XK = F2(-6, 6)


def shift(cs, c):
    out = [F2(0)] * len(cs)
    for j, cj in enumerate(cs):
        for i in range(j + 1):
            out[i] = out[i] + cj * comb(j, i) * (c ** (j - i))
    return out


def poly(cs, x):
    r = F2(0)
    for c in reversed(cs):
        r = r * x + c
    return r


assert shift(XDEN, XK) == [F2(0), F2(0), F2(1)], "x_den is not t^2"
assert shift(YDEN, XK) == [F2(0), F2(0), F2(0), F2(1)], "y_den is not t^3"
a0, a1, a2, a3 = shift(XNUM, XK)
b0, b1, b2, b3 = shift(YNUM, XK)
c = b3 / a3
assert a3.b == 0 and c.b == 0 and b2 == F2(0)
assert a2 / a3 == XK and a1 / a3 == F2(0, 48) and a0 / a3 == F2(16, 16)
assert b1 == -(c * a1) and b0 == -(c * a0 * 2)


def iso_rfc(x, y):
    return poly(XNUM, x) / poly(XDEN, x), y * poly(YNUM, x) / poly(YDEN, x)


def iso_velu(x, y):
    s = F2(4) / (x - XK)
    X = a3 * (x + F2(0, 12) * s + F2(1, 1) * s * s)
    F = c * a3 / 2 * (F2(2) - F2(0, 6) * s * s - F2(1, 1) * s * s * s)
    return X, y * F, F


random.seed(2026)
for _ in range(64):
    x = F2(random.randrange(p), random.randrange(p))
    y = F2(random.randrange(p), random.randrange(p))
    g = x * x * x + A * x + B
    X, Y, F = iso_velu(x, y)
    Xr, Yr = iso_rfc(x, y)
    assert X == Xr and Y == Yr, "velu form != RFC map"
    assert Y * Y - X * X * X - BE == F * F * (y * y - g), "curve identity"

C = -B / A
K = 253 * pow(60, p - 2, p) % p
assert C == F2(K) * F2(-1, 1)
xi = F2(-2, -1)
for _ in range(64):
    u = F2(random.randrange(p), random.randrange(p))
    tv1 = xi * u * u
    tv2 = tv1 * tv1 + tv1
    x1 = C * (F2(1) + tv2.inv())
    x2 = tv1 * x1
    assert (x1 - C) * tv2 == C
    assert (x2 * 60 - F2(253) * F2(-1, 1) * tv1) * (tv1 + 1) == F2(253) * F2(-1, 1)


def limbs(v):
    return ", ".join(f"0x{(v >> (64 * i)) & (2**64 - 1):016x}" for i in range(6))


print("velu form, curve identity and pins verified")
print(f"a3          = [{limbs(a3.a)}]")
print(f"c a3 / 2    = [{limbs((c * a3 / 2).a)}]")
print(f"K = 253/60  = [{limbs(K)}]")
print(f"x_k.c0      = [{limbs(XK.a)}]")
