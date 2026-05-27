from pairing_utils import toy_curve

G, Fp, e = toy_curve()

print("=== Basic ===")
print(f"G(3) = {G(3)}")
print(f"G(2) + G(3) = {G(2) + G(3)}")
print(f"5 * G(1) = {5 * G(1)}")

print("\n=== Bilinearity ===")
print(f"e(G(3), G(5)) == e(G(1), G(1))^15: {e(G(3), G(5)) == e(G(1), G(1))**15}")
print(f"e(G(3), G(5)) == e(G(1), G(15)): {e(G(3), G(5)) == e(G(1), G(15))}")

print("\n=== sum() ===")
total = sum([G(i) for i in range(1, 6)], G(0))
print(f"sum(G(1..5)) == G(15): {total == G(15)}")

print(Fp)
