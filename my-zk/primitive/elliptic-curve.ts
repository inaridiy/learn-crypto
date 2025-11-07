import { ExtendedFiniteFieldElement } from "./extended-finite-field";
import { FiniteField, FiniteFieldElement } from "./finite-field";
import { CyclicGroup, AbelianGroupElement, Field, FieldElement } from "./interface";
import { PolynomialOnFF } from "./polynomial";

export type ECPointLike<F extends FieldElement<F, FLike>, FLike = unknown> = Readonly<{
  x: F | FLike;
  y: F | FLike;
}>;

export class EllipticCurve<F extends FieldElement<F, FLike>, FLike = unknown> {
  public readonly field: Field<F, FLike>;

  public readonly doublingScalar: F;
  public readonly triplingScalar: F;

  constructor(public readonly a: F, public readonly b: F) {
    this.field = a.structure;

    const one = this.field.one();
    this.doublingScalar = one.scale(2n);
    this.triplingScalar = one.scale(3n);
  }

  isOnCurve(pointLike: ECPointLike<F, FLike>): boolean {
    const { x, y } =
      pointLike instanceof ECPoint
        ? pointLike
        : { x: this.field.from(pointLike.x), y: this.field.from(pointLike.y) };

    if (x.isZero() && y.isZero()) return true;
    const left = y.pow(2n);
    const right = x.pow(3n).add(this.a.mul(x)).add(this.b);
    return left.eq(right);
  }

  zero(): ECPoint<this, F, FLike> {
    return new ECPoint(this, {
      x: this.field.zero(),
      y: this.field.zero(),
    });
  }

  from(value: ECPoint<this, F, FLike> | ECPointLike<F, FLike>): ECPoint<this, F, FLike> {
    if (value instanceof ECPoint) return value.clone();

    const x = this.field.from(value.x);
    const y = this.field.from(value.y);
    return new ECPoint(this, { x, y });
  }
}

export class ECPointCyclicGroup<F extends FieldElement<F, FLike>, FLike = unknown>
  extends EllipticCurve<F, FLike>
  implements CyclicGroup<ECPoint<ECPointCyclicGroup<F, FLike>, F, FLike>, ECPointLike<F, FLike>>
{
  public readonly order: bigint;
  readonly #generator: ECPoint<ECPointCyclicGroup<F, FLike>, F, FLike>;

  constructor(
    public readonly a: F,
    public readonly b: F,
    generatorPoint: ECPointLike<F, FLike>,
    order: bigint
  ) {
    super(a, b);
    this.order = order;
    this.#generator = this.from(generatorPoint);
  }

  generator() {
    return this.#generator;
  }
}

export class ECPoint<
  S extends EllipticCurve<F, FLike>,
  F extends FieldElement<F, FLike>,
  FLike = unknown
> {
  public readonly structure: S;
  public readonly x: F;
  public readonly y: F;

  constructor(structure: S, point: { x: F; y: F }) {
    if (!structure.isOnCurve(point)) throw new Error("Point is not on curve");
    this.structure = structure;
    this.x = point.x;
    this.y = point.y;
  }

  clone(): ECPoint<S, F, FLike> {
    return new ECPoint(this.structure, { x: this.x, y: this.y });
  }

  eq(other: ECPoint<S, F, FLike>): boolean {
    return this.x.eq(other.x) && this.y.eq(other.y);
  }

  isZero(): boolean {
    return this.x.isZero() && this.y.isZero();
  }

  add(other: ECPoint<S, F, FLike>): ECPoint<S, F, FLike> {
    if (this.isZero()) return other;
    if (other.isZero()) return this;

    const { x: x1, y: y1 } = this;
    const { x: x2, y: y2 } = other;

    if (!x1.eq(x2)) {
      const slope = y1.sub(y2).div(x1.sub(x2));
      const x3 = slope.pow(2n).sub(x1).sub(x2);
      const y3 = slope.mul(x1.sub(x3)).sub(y1);
      return new ECPoint(this.structure, { x: x3, y: y3 });
    } else if (!y1.eq(y2)) {
      return this.structure.zero();
    } else {
      const slope = x1
        .pow(2n)
        .mul(this.structure.triplingScalar)
        .add(this.structure.a)
        .div(y1.mul(this.structure.doublingScalar));
      const x3 = slope.pow(2n).sub(x1.mul(this.structure.doublingScalar));
      const y3 = slope.mul(x1.sub(x3)).sub(y1);

      return new ECPoint(this.structure, { x: x3, y: y3 });
    }
  }

  sub(other: ECPoint<S, F, FLike>): ECPoint<S, F, FLike> {
    if (other.isZero()) return this;
    return this.add(other.negate());
  }

  negate(): ECPoint<S, F, FLike> {
    if (this.isZero()) return this.structure.zero();
    return new ECPoint(this.structure, {
      x: this.x,
      y: this.y.negate(),
    });
  }

  scale(n: bigint): ECPoint<S, F, FLike> {
    if (n < 0n) throw new Error("n must be positive");
    if (n === 0n) return this.structure.zero();

    let result = this.structure.zero();
    let point = this.clone();
    while (n > 0n) {
      if (n % 2n === 1n) result = point.add(result);
      point = point.add(point);
      n >>= 1n;
    }
    return result;
  }
}

export type ECPointOnFF = ECPoint<
  EllipticCurve<FiniteFieldElement, bigint>,
  FiniteFieldElement,
  bigint
>;

export type ECPointCyclicGroupOnFF = ECPoint<
  ECPointCyclicGroup<FiniteFieldElement, bigint>,
  FiniteFieldElement,
  bigint
>;

export type ECPointOnExtendedFF = ECPoint<
  EllipticCurve<ExtendedFiniteFieldElement, PolynomialOnFF>,
  ExtendedFiniteFieldElement,
  PolynomialOnFF
>;

export type ECPointCyclicGroupOnExtendedFF = ECPoint<
  ECPointCyclicGroup<ExtendedFiniteFieldElement, PolynomialOnFF>,
  ExtendedFiniteFieldElement,
  PolynomialOnFF
>;

if (import.meta.vitest) {
  const { describe, test, expect } = import.meta.vitest;
  const { FiniteField } = await import("./finite-field");

  // Use a small finite field for testing y^2 = x^3 + ax + b
  const p = 17n;
  const F = new FiniteField(p);

  // Curve y^2 = x^3 + 2x + 2 (mod 17)
  const a = F.from(2n);
  const b = F.from(2n);
  const curve = new EllipticCurve<FiniteFieldElement, bigint>(a, b); //dummy
  const O = curve.zero(); // Point at infinity (identity element)

  // Find some points on the curve y^2 = x^3 + 2x + 2 (mod 17)
  // x=0: y^2 = 2 (no solution)
  // x=1: y^2 = 1+2+2 = 5 (no solution)
  // x=2: y^2 = 8+4+2 = 14 (no solution)
  // x=3: y^2 = 27+6+2 = 35 = 1 (mod 17) => y = 1, 16
  // x=4: y^2 = 64+8+2 = 74 = 6 (mod 17) (no solution)
  // x=5: y^2 = 125+10+2 = 137 = 1 (mod 17) => y = 1, 16
  // x=6: y^2 = 216+12+2 = 230 = 9 (mod 17) => y = 3, 14

  const P1 = curve.from({ x: 3n, y: 1n });
  const P2 = curve.from({ x: 5n, y: 1n });
  const P3 = curve.from({ x: 6n, y: 3n });
  const P1_neg = curve.from({ x: 3n, y: 16n });

  describe("EllipticCurve", () => {
    test("isOnCurve", () => {
      expect(curve.isOnCurve({ x: F.from(3n), y: F.from(1n) })).toBe(true);
      expect(curve.isOnCurve({ x: F.from(3n), y: F.from(16n) })).toBe(true);
      expect(curve.isOnCurve({ x: F.from(5n), y: F.from(1n) })).toBe(true);
      expect(curve.isOnCurve({ x: F.from(6n), y: F.from(3n) })).toBe(true);
      expect(curve.isOnCurve({ x: F.from(1n), y: F.from(1n) })).toBe(false); // y^2=5 != 1^2=1
      expect(curve.isOnCurve(O)).toBe(true); // Point at infinity is on curve
    });

    test("factory methods", () => {
      expect(O.isZero()).toBe(true);
      const p = curve.from({ x: F.from(3n), y: F.from(1n) });
      expect(p.x.eq(F.from(3n))).toBe(true);
      expect(p.y.eq(F.from(1n))).toBe(true);
      expect(() => curve.from({ x: F.from(1n), y: F.from(1n) })).toThrow("Point is not on curve");
    });

    test("eq", () => {
      const P1_clone = curve.from({ x: F.from(3n), y: F.from(1n) });
      expect(P1.eq(P2)).toBe(false);
      expect(P1.eq(P1_clone)).toBe(true);
      expect(P1.eq(O)).toBe(false);
      expect(O.eq(O)).toBe(true);
    });

    test("neg", () => {
      expect(P1.negate().eq(P1_neg)).toBe(true);
      expect(P1_neg.negate().eq(P1)).toBe(true);
      expect(O.negate().eq(O)).toBe(true);
    });

    test("add (point addition)", () => {
      // P1 + O = P1
      expect(P1.add(O).eq(P1)).toBe(true);
      // O + P1 = P1
      expect(O.add(P1).eq(P1)).toBe(true);
      // P1 + (-P1) = O
      expect(P1.add(P1_neg).eq(O)).toBe(true);

      // P1(3,1) + P2(5,1) : Different x, Different y (but y1=y2 here)
      // m = (1-1)/(3-5) -> leads to division by zero if not handled for y1=y2
      // However, the formula handles y1=y2 case separately from point doubling.
      // Let's re-verify the case: If x1 != x2 and y1 == y2, m = 0.
      // m = (1-1)/(3-5) = 0 / -2 = 0
      // x3 = m^2 - x1 - x2 = 0^2 - 3 - 5 = -8 = 9 (mod 17)
      // y3 = m(x1 - x3) - y1 = 0(3 - 9) - 1 = -1 = 16 (mod 17)
      const expected_P1P2 = curve.from({ x: F.from(9n), y: F.from(16n) });
      expect(P1.add(P2).eq(expected_P1P2)).toBe(true);

      // P1(3,1) + P3(6,3) : Different x, Different y
      // m = (1-3)/(3-6) = -2 / -3 = 2 * inv(3)
      // inv(3) mod 17: 3*6 = 18 = 1. inv(3) = 6.
      // m = 2 * 6 = 12
      // x3 = m^2 - x1 - x3 = 12^2 - 3 - 6 = 144 - 9 = 135
      // 135 mod 17: 17 * 7 = 119, 135 - 119 = 16. x3 = 16
      // y3 = m(x1 - x3) - y1 = 12(3 - 16) - 1 = 12 * (-13) - 1 = 12 * 4 - 1
      // y3 = 48 - 1 = 47
      // 47 mod 17: 17 * 2 = 34, 47 - 34 = 13. y3 = 13
      const expected_P1P3 = curve.from({ x: F.from(16n), y: F.from(13n) });
      expect(P1.add(P3).eq(expected_P1P3)).toBe(true);
    });

    test("add (point doubling)", () => {
      // P1(3,1) + P1(3,1)
      // m = (3*x1^2 + a) / (2*y1)
      // m = (3*3^2 + 2) / (2*1) = (3*9 + 2) / 2 = (27 + 2) / 2 = 29 / 2
      // m = 12 / 2 (mod 17)
      // inv(2) mod 17: 2*9 = 18 = 1. inv(2) = 9.
      // m = 12 * 9 = 108
      // 108 mod 17: 17 * 6 = 102, 108 - 102 = 6. m = 6
      // x3 = m^2 - 2*x1 = 6^2 - 2*3 = 36 - 6 = 30 = 13 (mod 17)
      // y3 = m(x1 - x3) - y1 = 6(3 - 13) - 1 = 6*(-10) - 1 = 6*7 - 1
      // y3 = 42 - 1 = 41
      // 41 mod 17: 17*2 = 34, 41 - 34 = 7. y3 = 7
      const expected_2P1 = curve.from({ x: F.from(13n), y: F.from(7n) });
      expect(P1.add(P1).eq(expected_2P1)).toBe(true);
    });

    test("sub", () => {
      // P1 - P2 = P1 + (-P2)
      const P2_neg = P2.negate(); // (5, 16)
      const expected_sub = P1.add(P2_neg);
      expect(P1.sub(P2).eq(expected_sub)).toBe(true);
      // P1 - P1 = O
      expect(P1.sub(P1).eq(O)).toBe(true);
      // P1 - O = P1
      expect(P1.sub(O).eq(P1)).toBe(true);
    });

    test("mul (scalar multiplication)", () => {
      // 0 * P1 = O
      expect(P1.scale(0n).eq(O)).toBe(true);
      // 1 * P1 = P1
      expect(P1.scale(1n).eq(P1)).toBe(true);
      // 2 * P1 = P1 + P1
      const expected_2P1 = P1.add(P1);
      expect(P1.scale(2n).eq(expected_2P1)).toBe(true);
      // 3 * P1 = P1 + P1 + P1 = (2*P1) + P1
      const expected_3P1 = expected_2P1.add(P1);
      expect(P1.scale(3n).eq(expected_3P1)).toBe(true);

      // Test multiplication by a larger number
      // 4 * P1 = 2 * (2*P1)
      const P_2P1 = P1.scale(2n);
      const P_4P1_method1 = P_2P1.add(P_2P1);
      const P_4P1_method2 = P1.scale(4n);
      expect(P_4P1_method2.eq(P_4P1_method1)).toBe(true);

      // Test with negative n (should throw error)
      expect(() => P1.scale(-1n)).toThrow("n must be positive");
    });
  });
}
