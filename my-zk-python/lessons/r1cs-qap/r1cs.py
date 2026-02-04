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

Aw = A * w
Bw = B * w
Cw = C * w

lhs = vector(F, [Aw[i] * Bw[i] for i in range(Aw.length())])  # Hadamard product

print("\nA =\n", A)
print("B =\n", B)
print("C =\n", C)
print("\nA*w =", Aw)
print("B*w =", Bw)
print("C*w =", Cw)
print("(A*w) ⊙ (B*w) =", lhs)
print("\nSatisfied?", lhs == Cw)
