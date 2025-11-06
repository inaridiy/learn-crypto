import { FiniteField } from "./finite-field";
import {
  CommutativeRing,
  CommutativeRingElement,
  EuclideanRingElement,
  Field,
  FieldElement,
} from "./interface";
import { pow } from "./math";

export class PolynomialFactory<T extends FieldElement<T, TLike>, TLike>
  implements CommutativeRing<Polynomial<T, TLike>, T[] | TLike[]>
{
  constructor(public readonly coeffField: Field<T, TLike>) {}

  zero() {
    return new Polynomial(this, [this.coeffField.zero()]);
  }

  one() {
    return new Polynomial(this, [this.coeffField.one()]);
  }

  from(value: T[] | TLike[] | Polynomial<T, TLike>) {
    if (value instanceof Polynomial) return value;
    return new Polynomial(this, value.map(this.coeffField.from));
  }
}

export class Polynomial<T extends FieldElement<T, TLike>, TLike>
  implements
    CommutativeRingElement<Polynomial<T, TLike>, T[]>,
    EuclideanRingElement<Polynomial<T, TLike>, T[]>
{
  public readonly structure: PolynomialFactory<T, TLike>;

  public readonly coeffs: T[];

  constructor(factory: PolynomialFactory<T, TLike>, coeffs: T[]) {
    this.structure = factory;

    const normalized = [...coeffs];
    while (normalized.length > 0 && normalized[normalized.length - 1].isZero()) {
      normalized.pop();
    }
    if (normalized.length === 0) {
      this.coeffs = [this.structure.coeffField.zero()];
    } else {
      this.coeffs = normalized;
    }
  }

  degree(): number {
    if (this.isZero()) return -1; //これでいいっけ
    return this.coeffs.length - 1;
  }

  eval(x: T) {
    let result = this.structure.coeffField.zero();
    for (let i = this.coeffs.length - 1; i >= 0; i--) {
      result = result.mul(x).add(this.coeffs[i]);
    }
    return result;
  }

  eq(other: Polynomial<T, TLike>) {
    if (this.coeffs.length !== other.coeffs.length) return false;
    for (let i = 0; i < this.coeffs.length; i++) {
      if (!this.coeffs[i].eq(other.coeffs[i])) return false;
    }
    return true;
  }

  clone() {
    return this.structure.from(this.coeffs.map((c) => c.clone()));
  }

  isZero(): boolean {
    return this.coeffs.length === 1 && this.coeffs[0].isZero();
  }

  isOne(): boolean {
    return this.coeffs.length === 1 && this.coeffs[0].isOne();
  }

  add(other: Polynomial<T, TLike>) {
    const maxLength = Math.max(this.coeffs.length, other.coeffs.length);
    const result: T[] = [];

    for (let i = 0; i < maxLength; i++) {
      const a = this.coeffs[i] ?? this.structure.coeffField.zero();
      const b = other.coeffs[i] ?? this.structure.coeffField.zero();
      result.push(a.add(b));
    }

    return this.structure.from(result);
  }

  sub(other: Polynomial<T, TLike>) {
    const maxLength = Math.max(this.coeffs.length, other.coeffs.length);
    const result: T[] = [];

    for (let i = 0; i < maxLength; i++) {
      const a = this.coeffs[i] ?? this.structure.coeffField.zero();
      const b = other.coeffs[i] ?? this.structure.coeffField.zero();
      result.push(a.sub(b));
    }
    return this.structure.from(result);
  }

  mul(other: Polynomial<T, TLike>) {
    if (this.isZero() || other.isZero()) return this.structure.zero();

    const resultDegree = this.degree() + other.degree();
    const result = Array(resultDegree + 1).fill(this.structure.coeffField.zero());

    for (let i = 0; i < this.coeffs.length; i++) {
      if (this.coeffs[i].isZero()) continue;
      for (let j = 0; j < other.coeffs.length; j++) {
        result[i + j] = result[i + j].add(this.coeffs[i].mul(other.coeffs[j]));
      }
    }

    return this.structure.from(result);
  }

  scale(n: bigint) {
    return this.structure.from(this.coeffs.map((c) => c.scale(n)));
  }

  negate() {
    return this.structure.from(this.coeffs.map((c) => c.negate()));
  }

  pow(n: bigint) {
    return pow<Polynomial<T, TLike>>(this, n);
  }

  // 体上の多項式固有の演算？

  quotient(other: Polynomial<T, TLike>): Polynomial<T, TLike> {
    if (other.isZero()) throw new Error("Division by zero");
    if (this.isZero()) return this.structure.zero();
    if (this.degree() < other.degree()) return this.structure.zero();

    if (other.degree() === 0) {
      const inv = other.coeffs[0].inverse();
      const coeffes = this.coeffs.map((c) => c.mul(inv));
      return this.structure.from(coeffes);
    }

    let quotient = this.structure.zero();
    let remainder = this.clone();
    const divisorLeadCoeffInv = other.coeffs[other.degree()].inverse();

    while (!remainder.isZero() && remainder.degree() >= other.degree()) {
      const degreeDiff = remainder.degree() - other.degree();
      const coeff = remainder.coeffs[remainder.degree()].mul(divisorLeadCoeffInv);

      const monomialCoefficients = Array(degreeDiff + 1).fill(this.structure.coeffField.zero());
      monomialCoefficients[degreeDiff] = coeff;
      const monomial = this.structure.from(monomialCoefficients);

      quotient = quotient.add(monomial);
      remainder = remainder.sub(monomial.mul(other));
    }

    return quotient;
  }

  remainder(other: Polynomial<T, TLike>): Polynomial<T, TLike> {
    if (other.isZero()) throw new Error("Modulo by zero polynomial");
    if (this.isZero()) return this.structure.zero();
    if (this.degree() < other.degree()) return this.clone();

    if (other.degree() === 0) return this.structure.zero();

    let remainder = this.clone();
    const divisorLeadCoeffInv = other.coeffs[other.degree()].inverse();

    while (!remainder.isZero() && remainder.degree() >= other.degree()) {
      const degreeDiff = remainder.degree() - other.degree();
      const coeff = remainder.coeffs[remainder.degree()].mul(divisorLeadCoeffInv);
      const monomialCoefficients = Array(degreeDiff + 1).fill(this.structure.coeffField.zero());
      monomialCoefficients[degreeDiff] = coeff;
      const monomial = this.structure.from(monomialCoefficients);
      remainder = remainder.sub(monomial.mul(other));
    }

    return remainder;
  }

  divmod(other: Polynomial<T, TLike>) {
    const q = this.quotient(other);
    const r = this.remainder(other);
    return [q, r] as const;
  }

  div(other: Polynomial<T, TLike>) {
    return this.quotient(other);
  }

  mod(other: Polynomial<T, TLike>) {
    return this.remainder(other);
  }
}

if (import.meta.vitest) {
  const { describe, test, expect } = import.meta.vitest;
  const F = new FiniteField(13n);
  const PF = new PolynomialFactory(F);

  describe("Polynomial", () => {
    const p1 = PF.from([1n, 2n, 3n]); // 3x^2 + 2x + 1
    const p2 = PF.from([4n, 5n]); // 5x + 4
    const zero = PF.zero();
    const one = PF.one();

    test("factory and basic properties", () => {
      expect(zero.coeffs.length).toBe(1);
      expect(zero.coeffs[0].eq(F.zero())).toBe(true);
      expect(one.coeffs.length).toBe(1);
      expect(one.coeffs[0].eq(F.one())).toBe(true);
      expect(p1.degree()).toBe(2);
      expect(p2.degree()).toBe(1);
      expect(zero.degree()).toBe(-1);
      expect(PF.from([0n]).isZero()).toBe(true);
      expect(PF.from([1n, 0n, 0n]).coeffs.length).toBe(1); // Normalization
    });

    test("eq", () => {
      const p1_clone = PF.from([1n, 2n, 3n]);
      expect(p1.eq(p2)).toBe(false);
      expect(p1.eq(p1_clone)).toBe(true);
    });

    test("add", () => {
      // (3x^2 + 2x + 1) + (5x + 4) = 3x^2 + 7x + 5
      const expected = PF.from([5n, 7n, 3n]);
      expect(p1.add(p2).eq(expected)).toBe(true);
      expect(p1.add(zero).eq(p1)).toBe(true);
    });

    test("sub", () => {
      // (3x^2 + 2x + 1) - (5x + 4) = 3x^2 - 3x - 3 = 3x^2 + 10x + 10 (mod 13)
      const expected = PF.from([10n, 10n, 3n]);
      expect(p1.sub(p2).eq(expected)).toBe(true);
      expect(p1.sub(p1).eq(zero)).toBe(true);
    });

    test("mul", () => {
      // (3x^2 + 2x + 1) * (5x + 4)
      // = 15x^3 + 12x^2 + 10x^2 + 8x + 5x + 4
      // = 15x^3 + 22x^2 + 13x + 4
      // = 2x^3 + 9x^2 + 0x + 4 (mod 13)
      const expected = PF.from([4n, 0n, 9n, 2n]);
      expect(p1.mul(p2).eq(expected)).toBe(true);
      expect(p1.mul(one).eq(p1)).toBe(true);
      expect(p1.mul(zero).eq(zero)).toBe(true);
    });

    test("scale", () => {
      // (3x^2 + 2x + 1) * 3 = 9x^2 + 6x + 3
      const expected = PF.from([3n, 6n, 9n]);
      expect(p1.scale(3n).eq(expected)).toBe(true);
    });

    test("div", () => {
      const p_dividend = PF.from([4n, 0n, 9n, 2n]); // 2x^3 + 9x^2 + 4
      const p_divisor = PF.from([4n, 5n]); // 5x + 4
      const expected_quotient = PF.from([1n, 2n, 3n]); // 3x^2 + 2x + 1
      const quotient = p_dividend.div(p_divisor);
      // Due to potential floating point issues in manual calculation, verify by multiplication
      // quotient * p_divisor + remainder = p_dividend
      const remainder = p_dividend.mod(p_divisor);
      const check = quotient.mul(p_divisor).add(remainder);
      // Use a direct check for division result (quotient)
      const q = p_dividend.div(p_divisor);
      expect(q.eq(expected_quotient)).toBe(true); // Direct check

      // Test division by constant
      const p_const = PF.from([2n]); // Constant polynomial 2
      const expected_div_const = PF.from([7n, 1n, 8n]); // (3x^2 + 2x + 1) / 2 mod 13 = (3x^2+2x+1)*inv(2)= (3x^2+2x+1)*7 = 21x^2+14x+7 = 8x^2+x+7
      expect(p1.div(p_const).eq(expected_div_const)).toBe(true);

      expect(() => p1.div(zero)).toThrow("Division by zero");
      expect(p1.div(p1).eq(one)).toBe(true);
    });

    test("mod", () => {
      const p_dividend = PF.from([5n, 7n, 3n, 1n]); // x^3 + 3x^2 + 7x + 5
      const p_divisor = PF.from([1n, 1n]); // x + 1
      // Remainder should be P(-1)
      // P(-1) = (-1)^3 + 3(-1)^2 + 7(-1) + 5 = -1 + 3 - 7 + 5 = 0
      const expected_remainder = PF.zero();
      expect(p_dividend.mod(p_divisor).eq(expected_remainder)).toBe(true);

      const p3 = PF.from([1n, 0n, 1n]); // x^2 + 1
      const p4 = PF.from([1n, 1n]); // x + 1
      // (x^2+1) mod (x+1) => (-1)^2 + 1 = 1+1 = 2
      const expected_rem2 = PF.from([2n]);
      expect(p3.mod(p4).eq(expected_rem2)).toBe(true);

      expect(() => p1.mod(zero)).toThrow("Modulo by zero polynomial");
      expect(p1.mod(p1).eq(zero)).toBe(true);
    });

    test("pow", () => {
      // (x+1)^2 = x^2 + 2x + 1
      const px1 = PF.from([1n, 1n]);
      const expected_pow2 = PF.from([1n, 2n, 1n]);
      expect(px1.pow(2n).eq(expected_pow2)).toBe(true);
      expect(p1.pow(0n).eq(one)).toBe(true);
      expect(p1.pow(1n).eq(p1)).toBe(true);
    });

    test("eval", () => {
      // P(x) = 3x^2 + 2x + 1
      // P(2) = 3*(2^2) + 2*2 + 1 = 3*4 + 4 + 1 = 12 + 4 + 1 = 17 = 4 mod 13
      const x_val = F.from(2n);
      const expected_eval = F.from(4n);
      expect(p1.eval(x_val).eq(expected_eval)).toBe(true);
      // P(0) = 1
      expect(p1.eval(F.zero()).eq(F.from(1n))).toBe(true);
    });
  });
}
