import { FiniteFieldElement } from "./finite-field";
import { Field, FieldElement } from "./interface";
import { extendedGCD, pow } from "./math";
import { Polynomial } from "./polynomial";

type PolyOnFF = Polynomial<FiniteFieldElement, bigint | FiniteFieldElement>;

export class ExtendedFiniteField implements Field<ExtendedFiniteFieldElement, PolyOnFF> {
  constructor(public readonly p: PolyOnFF) {}

  zero() {
    return new ExtendedFiniteFieldElement(this, this.p, this.p.structure.zero());
  }

  one() {
    return new ExtendedFiniteFieldElement(this, this.p, this.p.structure.one());
  }

  from(value: ExtendedFiniteFieldElement | PolyOnFF) {
    if (value instanceof ExtendedFiniteFieldElement) return value;
    return new ExtendedFiniteFieldElement(this, this.p, value);
  }
}

export class ExtendedFiniteFieldElement
  implements FieldElement<ExtendedFiniteFieldElement, PolyOnFF>
{
  public readonly structure: ExtendedFiniteField;

  public readonly CommutativeRingElement = true;
  public readonly FieldElement = true;

  public readonly p: PolyOnFF;
  public readonly n: PolyOnFF;

  constructor(structure: ExtendedFiniteField, p: PolyOnFF, n: PolyOnFF) {
    this.structure = structure;
    this.p = p;
    this.n = n.mod(this.p);
  }

  eq(other: ExtendedFiniteFieldElement): boolean {
    return this.n.eq(other.n);
  }

  clone(): ExtendedFiniteFieldElement {
    return this.structure.from(this.n.clone());
  }

  isZero(): boolean {
    return this.n.isZero();
  }

  isOne(): boolean {
    return this.n.isOne();
  }

  add(other: ExtendedFiniteFieldElement) {
    return this.structure.from(this.n.add(other.n));
  }

  sub(other: ExtendedFiniteFieldElement) {
    return this.structure.from(this.n.sub(other.n));
  }

  mul(other: ExtendedFiniteFieldElement) {
    return this.structure.from(this.n.mul(other.n));
  }

  scale(n: bigint) {
    return this.structure.from(this.n.scale(n));
  }

  inverse() {
    const [g, x, _] = extendedGCD(this.n, this.p);
    if (g.degree() !== 0) throw new Error("Not invertible");
    else if (g.coeffs[0].n !== 1n) return this.structure.from(x.div(g));
    else return this.structure.from(x);
  }

  div(other: ExtendedFiniteFieldElement) {
    return this.mul(other.inverse());
  }

  pow(n: bigint) {
    return pow<ExtendedFiniteFieldElement>(this, n);
  }
}

if (import.meta.vitest) {
  const { describe, test, expect } = import.meta.vitest;
  const { FiniteField } = await import("./finite-field");
  const { PolynomialFactory } = await import("./polynomial");

  // GF(2^3) over GF(2) using irreducible polynomial x^3 + x + 1
  const F2 = new FiniteField(2n);
  const PF2 = new PolynomialFactory(F2);
  const modPoly = PF2.from([1n, 1n, 0n, 1n]); // x^3 + x + 1
  const EFF = new ExtendedFiniteField(modPoly);

  // Elements represented as polynomials in 'a', where a^3 + a + 1 = 0
  const zero = EFF.zero(); // 0
  const one = EFF.one(); // 1
  const a = EFF.from(PF2.from([0n, 1n])); // a
  const a2 = EFF.from(PF2.from([0n, 0n, 1n])); // a^2
  const a_plus_1 = EFF.from(PF2.from([1n, 1n])); // a + 1
  const a2_plus_a = EFF.from(PF2.from([0n, 1n, 1n])); // a^2 + a
  const a2_plus_a_plus_1 = EFF.from(PF2.from([1n, 1n, 1n])); // a^2 + a + 1
  const a2_plus_1 = EFF.from(PF2.from([1n, 0n, 1n])); // a^2 + 1

  describe("ExtendedFiniteField GF(2^3)", () => {
    test("factory methods", () => {
      expect(zero.n.isZero()).toBe(true);
      expect(one.n.eq(PF2.one())).toBe(true);
      // Test 'from' and normalization (mod)
      // x^3 = x + 1 (mod x^3+x+1 over F2)
      const a3 = EFF.from(PF2.from([0n, 0n, 0n, 1n]));
      expect(a3.eq(a_plus_1)).toBe(true);
    });

    test("eq", () => {
      const a_clone = EFF.from(PF2.from([0n, 1n]));
      expect(a.eq(a2)).toBe(false);
      expect(a.eq(a_clone)).toBe(true);
    });

    test("isZero", () => {
      expect(zero.isZero()).toBe(true);
      expect(one.isZero()).toBe(false);
      expect(a.isZero()).toBe(false);
    });

    test("add", () => {
      // a + a^2
      expect(a.add(a2).eq(a2_plus_a)).toBe(true);
      // (a+1) + (a^2+a) = a^2 + 2a + 1 = a^2 + 1 (mod 2)
      expect(a_plus_1.add(a2_plus_a).eq(a2_plus_1)).toBe(true);
      // a + a = 2a = 0 (mod 2)
      expect(a.add(a).eq(zero)).toBe(true);
      expect(a.add(zero).eq(a)).toBe(true);
    });

    test("sub", () => {
      // Subtraction is the same as addition in GF(2^n)
      expect(a.sub(a2).eq(a2_plus_a)).toBe(true);
      expect(a_plus_1.sub(a2_plus_a).eq(a2_plus_1)).toBe(true);
      expect(a.sub(a).eq(zero)).toBe(true);
      expect(a.sub(zero).eq(a)).toBe(true);
    });

    test("mul", () => {
      // a * a = a^2
      expect(a.mul(a).eq(a2)).toBe(true);
      // a * a^2 = a^3 = a + 1 (mod x^3+x+1)
      expect(a.mul(a2).eq(a_plus_1)).toBe(true);
      // (a+1) * (a^2+1) = a^3 + a + a^2 + 1 = (a+1) + a + a^2 + 1 = a^2 (mod 2)
      expect(a_plus_1.mul(a2_plus_1).eq(a2)).toBe(true);
      expect(a.mul(one).eq(a)).toBe(true);
      expect(a.mul(zero).eq(zero)).toBe(true);
    });

    test("scale (field is F2^n, so scaling by bigint n means repeated addition)", () => {
      // Scale by 0 -> zero
      expect(a.scale(0n).eq(zero)).toBe(true);
      // Scale by 1 -> itself
      expect(a.scale(1n).eq(a)).toBe(true);
      // Scale by 2 (even) -> zero (since 1+1=0 in F2)
      expect(a.scale(2n).eq(zero)).toBe(true);
      // Scale by 3 (odd) -> itself (since 1+1+1=1 in F2)
      expect(a.scale(3n).eq(a)).toBe(true);
    });

    test("pow", () => {
      // a^2
      expect(a.pow(2n).eq(a2)).toBe(true);
      // a^3 = a+1
      expect(a.pow(3n).eq(a_plus_1)).toBe(true);
      // a^0 = 1
      expect(a.pow(0n).eq(one)).toBe(true);
      // (a+1)^2 = a^2 + 2a + 1 = a^2 + 1
      expect(a_plus_1.pow(2n).eq(a2_plus_1)).toBe(true);
      // Field is GF(2^3), so order is 2^3-1=7. a^7 should be 1.
      expect(a.pow(7n).eq(one)).toBe(true);
    });
  });
}
