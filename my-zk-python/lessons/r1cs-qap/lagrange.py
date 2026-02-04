from sage.all import *

R = PolynomialRing(QQ, "x")
x = R.gen()
p = R.lagrange_polynomial([(QQ(1), QQ(5)), (QQ(2), QQ(0)), (QQ(3), QQ(3))])

print(p)
