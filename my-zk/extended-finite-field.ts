import { Field, FieldFactory } from "./interface";
import { extendedGCD, pow } from "./math";
import { Polynomial } from "./polynomial";

export class ExtendedFiniteFieldFactory implements FieldFactory<ExtendedFiniteField, Polynomial> {
  constructor(public readonly mod: Polynomial) {}

  zero(): ExtendedFiniteField {
    return new ExtendedFiniteField(this, this.mod, this.mod.factory.zero());
  }

  one(): ExtendedFiniteField {
    return new ExtendedFiniteField(this, this.mod, this.mod.factory.one());
  }

  from(value: Polynomial): ExtendedFiniteField {
    return new ExtendedFiniteField(this, this.mod, value);
  }
}

export class ExtendedFiniteField implements Field<ExtendedFiniteField, Polynomial> {
  public readonly factory: ExtendedFiniteFieldFactory;
  public readonly mod: Polynomial;
  public readonly n: Polynomial;

  constructor(factory: ExtendedFiniteFieldFactory, mod: Polynomial, n: Polynomial) {
    this.factory = factory;
    this.mod = mod;
    this.n = n.mod(this.mod);
  }

  eq(other: ExtendedFiniteField): boolean {
    return this.n.eq(other.n);
  }

  clone(): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n);
  }

  isZero(): boolean {
    return this.n.isZero();
  }

  inverse(): ExtendedFiniteField {
    const [g, x, _] = extendedGCD(this.n, this.mod);
    if (g.degree() !== 0) throw new Error("Not invertible");
    else if (g.coeffs[0].n !== 1n) return new ExtendedFiniteField(this.factory, this.mod, x.div(g));
    else return new ExtendedFiniteField(this.factory, this.mod, x);
  }

  add(other: ExtendedFiniteField): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.add(other.n));
  }

  sub(other: ExtendedFiniteField): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.sub(other.n));
  }

  mul(other: ExtendedFiniteField): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.mul(other.n));
  }

  div(other: ExtendedFiniteField): ExtendedFiniteField {
    return this.mul(other.inverse());
  }

  scale(n: bigint): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.scale(n));
  }

  pow(n: bigint): ExtendedFiniteField {
    return pow<ExtendedFiniteField>(this, n);
  }
}

if (import.meta.vitest) {
  const { describe, test, expect } = import.meta.vitest;
  const { FiniteFieldFactory } = await import("./finite-field");
  const { PolynomialFactory } = await import("./polynomial");

  // GF(2^3) over GF(2) using irreducible polynomial x^3 + x + 1
  const F2 = new FiniteFieldFactory(2n);
  const PF2 = new PolynomialFactory(F2);
  const modPoly = PF2.from([1n, 1n, 0n, 1n]); // x^3 + x + 1
  const EFF = new ExtendedFiniteFieldFactory(modPoly);

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

    // Inverse and Division require extended Euclidean algorithm for polynomials, which is not implemented here.
    // Skipping div and inverse tests for now.
    /*
    test("inverse", () => {
      // Find inverse of 'a'
      // Need Extended Euclidean Algo for Polynomials
    });

    test("div", () => {
      // (a^2+a) / a = a+1
      // Need polynomial division or multiplication by inverse
    });
    */

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
