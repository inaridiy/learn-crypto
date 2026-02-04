from sage.all import *


F = QQ

# ---- witness (x=2) ----
x = F(2)
i1 = F(5) * x * x
i2 = i1 * x
out = i2 + F(3) * x

w = vector(F, [F(1), x, out, i1, i2])
print("w = [1, x, out, i1, i2] =", w)

A = matrix(F, [[0, 5, 0, 0, 0], [0, 0, 0, 1, 0], [0, 3, 0, 0, 1]])
B = matrix(F, [[0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [1, 0, 0, 0, 0]])
C = matrix(F, [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0]])

R = PolynomialRing(F, "x")
t = R.gen()
m = A.nrows()  # 制約数 = 3
n = A.ncols()  # witness要素数 = 5

# 評価点
eval_points = [F(i + 1) for i in range(m)]  # [1, 2, 3]


def plot_poly(f, xmin, xmax, filename):
    g = f.plot(xmin=xmin, xmax=xmax)
    g.save(filename)
    return filename


# 各列をラグランジュ補間して多項式を得る
def matrix_to_polynomials(M):
    polys = []
    for j in range(n):
        points = [(eval_points[i], M[i, j]) for i in range(m)]
        polys.append(R.lagrange_polynomial(points))
    return vector(R, polys)


AP = matrix_to_polynomials(A)
BP = matrix_to_polynomials(B)
CP = matrix_to_polynomials(C)

print("AP =", AP)
print("BP =", BP)
print("CP =", CP)

AX = AP * w
BX = BP * w
CX = CP * w

print("AX =", AX)
print("BX =", BX)
print("CX =", CX)

PX = AX * BX - CX

print("P(x) =", PX)
print("P(1) =", PX(1))
print("P(2) =", PX(2))
print("P(3) =", PX(3))

plot_poly(PX, -1, 4, "/home/sage/project/my-zk-python/lessons/r1cs-qap/px.png")

tx = (t - 1)(t - 2)(t - 3)
print("t(x)", t)

hx = PX / tx
print("h(x)", hx)

print("isOK", hx * tx == PX)
