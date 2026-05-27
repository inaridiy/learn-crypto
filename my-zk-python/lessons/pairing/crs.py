from sage.all import *
from pairing_utils import toy_curve
from random import random
from math import floor

G, F, e = toy_curve()

R = PolynomialRing(F, "x")
x = R.gen()


# 次元数nを取って、[g,g^s,gs^{s^2}...]を返す
def generate_crs(n):
    s = F(floor(random() * F.order()))
    g = G(1)
    points = []
    for n in range(n):
        points.append(g)
        g = g * s

    return points, s


# 多項式とcrs([g,g^s,gs^{s^2}...])を取って、g^{f(s)}を返す
def eval_with_crs(fx, crs):
    coeffs = fx.list()
    acc = G(0)  # 単位元
    for a, g_si in zip(coeffs, crs):
        acc = acc + (g_si * a)
    return acc


f = x**4 + 2 * x**2 + 10

crs, s = generate_crs(10)

fs1 = eval_with_crs(f, crs)
fs2 = G(f(s))

print(s, f(s), F.order())
print("fs1 = ", fs1)
print("fs2 = ", fs2)
print("is ok", fs1 == fs2)
