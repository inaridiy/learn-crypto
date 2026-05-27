# pairing_utils.py
# Pinocchio教育用ペアリングユーティリティ（薄いラッパー）
# 実行: sage -python your_script.py

from sage.all import *


def setup_curve(p, A, B, n, k, tr, gen1_coords, gen2_coords, modulus_coeffs):
    """曲線セットアップ。(g1, g2, G, pairing) を返す"""
    n, k, tr = Integer(n), Integer(k), Integer(tr)
    Fp = GF(p)
    E = EllipticCurve(Fp, [A, B])

    R = Fp["x"]
    x = R.gen()
    modulus = sum(Integer(c) * x**i for i, c in enumerate(modulus_coeffs))
    Fpk = GF((p, k), modulus=modulus, names=("a",))
    a = Fpk.gen()

    EK = E.base_extend(Fpk)

    # G1生成元
    g1 = EK(Integer(gen1_coords[0]), Integer(gen1_coords[1]))

    # G2生成元
    x_coeffs, y_coeffs = gen2_coords
    g2_x = sum(Integer(c) * a**i for i, c in enumerate(x_coeffs))
    g2_y = sum(Integer(c) * a**i for i, c in enumerate(y_coeffs))
    g2 = EK(g2_x, g2_y)

    # Gファクトリ
    def G(scalar):
        return GPoint(Integer(scalar) * g1, Integer(scalar) * g2)

    # ペアリング関数
    def pairing(a, b):
        if isinstance(a, GPoint) and isinstance(b, GPoint):
            return a.g1.ate_pairing(b.g2, n, k, tr)
        elif isinstance(a, GPoint):
            return a.g1.ate_pairing(b, n, k, tr)
        elif isinstance(b, GPoint):
            return a.ate_pairing(b.g2, n, k, tr)
        else:
            return a.ate_pairing(b, n, k, tr)

    F = GF(g1.order())

    return G, F, pairing


class GPoint:
    """G1とG2を同時に保持"""

    def __init__(self, g1, g2):
        self.g1 = g1
        self.g2 = g2

    def __add__(self, other):
        return GPoint(self.g1 + other.g1, self.g2 + other.g2)

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __sub__(self, other):
        return GPoint(self.g1 - other.g1, self.g2 - other.g2)

    def __neg__(self):
        return GPoint(-self.g1, -self.g2)

    def __mul__(self, scalar):
        s = Integer(scalar)
        return GPoint(s * self.g1, s * self.g2)

    def __rmul__(self, scalar):
        return self * scalar

    def __eq__(self, other):
        return self.g1 == other.g1 and self.g2 == other.g2

    def __repr__(self):
        return f"G(g1={self.g1}, g2={self.g2})"


# プリセット: toy曲線
# --- toy_curve をこれに差し替え ---

_TOY_CACHE = None


def toy_curve():
    """
    BLS12-381 を使う “非toy” 曲線。
    sage -python で動くように、preparser 記法（R.<x>, ^）は使わない。
    戻り値: (G, F, pairing)  ※既存インターフェース維持
    """
    global _TOY_CACHE
    if _TOY_CACHE is not None:
        return _TOY_CACHE

    from sage.all import GF, EllipticCurve, Integer, PolynomialRing

    # BLS12-381 prime p と subgroup order r
    p = Integer(
        "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab"
    )
    r = Integer("0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001")

    # --- Field tower: Fp -> Fp2 = Fp[u]/(u^2+1) -> Fp12 = Fp2[w]/(w^6-(u+1)) ---
    Fp = GF(p)

    Rp = PolynomialRing(Fp, "X")
    X = Rp.gen()
    Fp2 = GF((p, 2), modulus=(X**2 + 1), names=("u",))
    u = Fp2.gen()

    Rp2 = PolynomialRing(Fp2, "Y")
    Y = Rp2.gen()
    Fp12 = Fp2.extension(Y**6 - (u + 1), names=("w",))
    w = Fp12.gen()

    # E: y^2 = x^3 + 4
    E = EllipticCurve(Fp12, [0, 4])

    # Twist side (G2): y^2 = x^3 + 4*(u+1) という形がよく使われる :contentReference[oaicite:1]{index=1}
    Et = EllipticCurve(Fp12, [0, 4 * (1 + u)])

    # --- Standard G1 generator (hex constants) ---
    g1x = Integer(
        "0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb"
    )
    g1y = Integer(
        "0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1"
    )
    g1 = E(Fp12(g1x), Fp12(g1y))

    # --- Standard G2 generator on Fp2 (hex constants; Sage Q&A に同じ値が載ってる) ---
    # x2 = x0 + x1*u, y2 = y0 + y1*u :contentReference[oaicite:2]{index=2}
    x0 = Integer(
        "0x024AA2B2F08F0A91260805272DC51051C6E47AD4FA403B02B4510B647AE3D1770BAC0326A805BBEFD48056C8C121BDB8"
    )
    x1 = Integer(
        "0x13E02B6052719F607DACD3A088274F65596BD0D09920B61AB5DA61BBDC7F5049334CF11213945D57E5AC7D055D042B7E"
    )
    y0 = Integer(
        "0x0CE5D527727D6E118CC9CDC6DA2E351AADFD9BAA8CBDD3A76D429A695160D12C923AC9CC3BACA289E193548608B82801"
    )
    y1 = Integer(
        "0x0606C4A02EA734CC32ACD2B02BC28B99CB3E287E85A763AF267492AB572E99AB3F370D275CEC1DA1AAA9075FF05F79BE"
    )

    x2 = Fp2(x0) + Fp2(x1) * u
    y2 = Fp2(y0) + Fp2(y1) * u

    # Twist 点として解釈（Fp12 へは自然に埋め込める）
    g2_tw = Et(Fp12(x2), Fp12(y2))

    # Sextic untwist: w^6 = (u+1) なので (x,y) -> (x/w^2, y/w^3) で E に写る
    g2 = E(Fp12(x2) / (w**2), Fp12(y2) / (w**3))

    # スカラー倍で (g1,g2) を同時に持つ
    def G(scalar):
        s = Integer(scalar)
        return GPoint(s * g1, s * g2)

    # スカラー体（Z_r）
    F = GF(r)

    def pairing(a, b):
        # 既存の便利インターフェース維持
        if isinstance(a, GPoint):
            a = a.g1
        if isinstance(b, GPoint):
            b = b.g2

        # ate_pairing は “Miller + naive exponentiation” で効率保証なし :contentReference[oaicite:3]{index=3}
        # Weil pairing は algorithm='pari' が選べる（PARI の ellweilpairing）:contentReference[oaicite:4]{index=4}
        try:
            return a.weil_pairing(b, r, algorithm="pari")
        except (TypeError, ValueError):
            # 古い Sage などで 'pari' が無い時のフォールバック
            return a.weil_pairing(b, r, algorithm="sage")

    _TOY_CACHE = (G, F, pairing)
    return _TOY_CACHE


if __name__ == "__main__":
    G, e = toy_curve()

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
