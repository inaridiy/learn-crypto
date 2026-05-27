from sage.all import *
from pairing_utils import toy_curve
from random import random
from math import floor

G, F, e = toy_curve()

R = PolynomialRing(F, "x")
t = R.gen()

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


x = F(2)
i1 = F(5) * x * x
i2 = i1 * x
out = i2 + F(3) * x

w = vector(F, [F(1), x, out, i1, i2])

A = matrix(F, [[0, 5, 0, 0, 0], [0, 0, 0, 1, 0], [0, 3, 0, 0, 1]])
B = matrix(F, [[0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [1, 0, 0, 0, 0]])
C = matrix(F, [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0]])

m = A.nrows()  # 制約数 = 3
n = A.ncols()  # witness要素数 = 5
eval_points = [F(i + 1) for i in range(m)]  # [1, 2, 3]

def matrix_to_polynomials(M):
    polys = []
    for j in range(n):
        points = [(eval_points[i], M[i, j]) for i in range(m)]
        polys.append(R.lagrange_polynomial(points))
    return vector(R, polys)



AP = matrix_to_polynomials(A)
BP = matrix_to_polynomials(B)
CP = matrix_to_polynomials(C)

# TTPの処理
crs,s = generate_crs(1000)
tx = (t - 1) * (t - 2) * (t - 3)
gts = eval_with_crs(tx,crs)

# 証明者の処理
AX = AP * w
BX = BP * w
CX = CP * w

PX = AX * BX - CX
hx = PX // tx

mid_start = 3

wio = w[:mid_start]
wmid = w[mid_start:]

gAmidS =  # g^{A_{mid}(s)}
gBmidS = # g^{B_{mid}(s)}
gCmidS = # g^{C_{mid}(s)}

ghs = # g^{h(s)}

print("検証者に wio,gAmidS,gBmidS,gBmidS,ghsを送る")
print("wio =",wio,wmid)
print("gAmidS =",gAmidS)
print("gBmidS =",gBmidS)
print("gCmidS =",gCmidS)
print("ghs =", ghs)

# 検証者の検証処理
APio = AP[:mid_start]
BPio = BP[:mid_start]
CPio = CP[:mid_start]

gAioS = # g^{A_{io}(s)}
gBioS = # g^{B_{io}(s)}
gCioS = # g^{C_{io}(s)}

gAs = # g^{A(s)}
gBs = # g^{B(s)}
gCs = # g^{C(s)}

lh = e(, ) # en(g^{A(s)},g^{B(s)})
# ペアリングの合成は加算ではなく乗算
rh = e(, ) * e(,G(1)) # en(g^{h(s)},g^{t(s)}) * en(g^{C(s)},g)

print("isOk",lh==rh)
