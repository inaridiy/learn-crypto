import { Field, FieldElement } from "./interface";
import { extendedGCDBigint, modBigint, pow } from "./math";

export class FiniteField implements Field<FiniteFieldElement, bigint> {
  constructor(public readonly p: bigint) {}

  zero() {
    return new FiniteFieldElement(this, this.p, 0n);
  }

  one() {
    return new FiniteFieldElement(this, this.p, 1n);
  }

  from(value: bigint | FiniteFieldElement) {
    if (value instanceof FiniteFieldElement) return value.clone();
    return new FiniteFieldElement(this, this.p, value);
  }
}

export class FiniteFieldElement implements FieldElement<FiniteFieldElement, bigint> {
  public readonly structure: FiniteField;

  public readonly p: bigint;
  public readonly n: bigint;

  constructor(structure: FiniteField, p: bigint, n: bigint) {
    this.structure = structure;
    this.p = p;
    this.n = modBigint(n, this.p);
  }

  eq(other: FiniteFieldElement): boolean {
    return this.n === other.n;
  }

  clone() {
    return new FiniteFieldElement(this.structure, this.p, this.n);
  }

  isZero(): boolean {
    return this.n === 0n;
  }

  isOne(): boolean {
    return this.n === 1n;
  }

  add(other: FiniteFieldElement) {
    return new FiniteFieldElement(this.structure, this.p, this.n + other.n);
  }

  sub(other: FiniteFieldElement) {
    return new FiniteFieldElement(this.structure, this.p, this.n - other.n);
  }

  mul(other: FiniteFieldElement) {
    return new FiniteFieldElement(this.structure, this.p, this.n * other.n);
  }

  scale(n: bigint) {
    return new FiniteFieldElement(this.structure, this.p, this.n * n);
  }

  negate() {
    return new FiniteFieldElement(this.structure, this.p, -this.n);
  }

  inverse() {
    const [gcd, x] = extendedGCDBigint(this.n, this.p);
    if (gcd !== 1n) throw new Error("No inverse");
    return new FiniteFieldElement(this.structure, this.p, x);
  }

  div(other: FiniteFieldElement) {
    if (other.n === 0n) throw new Error("Division by zero");
    const otherInv = other.inverse();
    return this.mul(otherInv);
  }

  quotient(other: FiniteFieldElement) {
    return this.div(other);
  }

  remainder(_: FiniteFieldElement) {
    return this.structure.zero();
  }

  divmod(other: FiniteFieldElement) {
    return [this.div(other), this.structure.zero()] as const;
  }

  pow(n: bigint) {
    return pow<FiniteFieldElement>(this, n);
  }
}

if (import.meta.vitest) {
  const { describe, test, expect } = import.meta.vitest;

  describe("FiniteField", () => {
    const p = 13n;
    const F = new FiniteField(p);
    const zero = F.zero();
    const one = F.one();
    const a = F.from(7n);
    const b = F.from(8n);
    const c = F.from(7n); // Same value as a

    test("factory methods", () => {
      expect(zero.n).toBe(0n);
      expect(one.n).toBe(1n);
      expect(F.from(15n).n).toBe(2n); // 15 mod 13 = 2
    });

    test("eq", () => {
      expect(a.eq(b)).toBe(false);
      expect(a.eq(c)).toBe(true);
      expect(a.eq(F.from(7n))).toBe(true);
    });

    test("clone", () => {
      const aClone = a.clone();
      expect(aClone.eq(a)).toBe(true);
      expect(aClone).not.toBe(a); // Ensure it's a different instance
    });

    test("isZero", () => {
      expect(zero.isZero()).toBe(true);
      expect(one.isZero()).toBe(false);
      expect(a.isZero()).toBe(false);
    });

    test("add", () => {
      const expected = F.from(2n); // (7 + 8) mod 13 = 15 mod 13 = 2
      expect(a.add(b).eq(expected)).toBe(true);
      expect(a.add(zero).eq(a)).toBe(true);
    });

    test("sub", () => {
      const expected = F.from(12n); // (7 - 8) mod 13 = -1 mod 13 = 12
      expect(a.sub(b).eq(expected)).toBe(true);
      const expected2 = F.from(1n); // (8 - 7) mod 13 = 1 mod 13 = 1
      expect(b.sub(a).eq(expected2)).toBe(true);
      expect(a.sub(zero).eq(a)).toBe(true);
      expect(a.sub(a).eq(zero)).toBe(true);
    });

    test("mul", () => {
      const expected = F.from(4n); // (7 * 8) mod 13 = 56 mod 13 = 4
      expect(a.mul(b).eq(expected)).toBe(true);
      expect(a.mul(one).eq(a)).toBe(true);
      expect(a.mul(zero).eq(zero)).toBe(true);
    });

    test("scale", () => {
      const expected = F.from(8n); // (7 * 3) mod 13 = 21 mod 13 = 8
      expect(a.scale(3n).eq(expected)).toBe(true);
      expect(a.scale(0n).eq(zero)).toBe(true);
      expect(a.scale(1n).eq(a)).toBe(true);
    });

    test("inverse", () => {
      const aInv = a.inverse(); // Inverse of 7 mod 13
      const expectedInv = F.from(2n); // 7 * 2 = 14 = 1 mod 13
      expect(aInv.eq(expectedInv)).toBe(true);
      expect(a.mul(aInv).eq(one)).toBe(true);
      expect(() => zero.inverse()).toThrow("No inverse"); // Inverse of 0 doesn't exist if gcd is not 1
      // Test inverse of 1
      expect(one.inverse().eq(one)).toBe(true);
    });

    test("div", () => {
      const expected = F.from(9n); // 7 / 8 = 7 * inv(8) = 7 * 5 = 35 = 9 mod 13. Let's recheck inv(8). 8 * 5 = 40 = 1 mod 13. Correct. Now 7 * 5 = 35 = 9 mod 13. Mistake in expected value.
      const bInv = b.inverse(); // inv(8) = 5
      expect(a.div(b).eq(expected)).toBe(true);
      expect(a.div(one).eq(a)).toBe(true);
      expect(() => a.div(zero)).toThrow("Division by zero");
      expect(a.div(a).eq(one)).toBe(true);
    });

    test("pow", () => {
      const expected = F.from(5n); // 7^3 mod 13 = 343 mod 13. 343 = 26 * 13 + 5. So 5.
      expect(a.pow(3n).eq(expected)).toBe(true);
      const expected2 = F.from(1n); // 7^0 mod 13 = 1
      expect(a.pow(0n).eq(one)).toBe(true); // pow(this, 0n) should return one()
      const expected3 = F.from(7n); // 7^1 mod 13 = 7
      expect(a.pow(1n).eq(a)).toBe(true);
      const expected4 = F.from(1n); // 7^12 mod 13 = 1 (Fermat's Little Theorem)
      expect(a.pow(12n).eq(one)).toBe(true);
      const expected5 = F.from(1n); // 0^0 mod 13 should be 1
      expect(zero.pow(0n).eq(one)).toBe(true);
      const expected6 = F.from(0n); // 0^3 mod 13 should be 0
      expect(zero.pow(3n).eq(zero)).toBe(true);
    });
  });
}
