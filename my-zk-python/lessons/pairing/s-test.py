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


hx = x**20 + 20 * x**9 + 23
tx = 9 * x**11 + 102 * x**2 + 10

px = hx * tx

crs, s = generate_crs(100)

gts = eval_with_crs(tx, crs)
ghs = eval_with_crs(hx, crs)
gps = eval_with_crs(px, crs)

lh = e(gts, ghs)
rh = e(gps, G(1))

print("isOk", lh == rh)
