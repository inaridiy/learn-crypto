import { FiniteField, FiniteFieldFactory } from "./finite-field";
import { Field, FieldFactory } from "./interface";
import { pow } from "./math";

export class PolynomialFactory implements FieldFactory<Polynomial, FiniteField[] | bigint[]> {
  constructor(public readonly coeffFactory: FiniteFieldFactory) {}

  zero(): Polynomial {
    return new Polynomial(this, [this.coeffFactory.zero()]);
  }

  one(): Polynomial {
    return new Polynomial(this, [this.coeffFactory.one()]);
  }

  from(coeffs: FiniteField[] | bigint[]): Polynomial {
    if (typeof coeffs[0] === "bigint") {
      return new Polynomial(
        this,
        coeffs.map((c) => this.coeffFactory.from(c as bigint))
      );
    }
    return new Polynomial(this, coeffs as FiniteField[]);
  }
}

export class Polynomial implements Field<Polynomial, FiniteField[]> {
  public readonly factory: PolynomialFactory;
  public readonly coeffFactory: FiniteFieldFactory;
  public readonly coeffs: FiniteField[];

  constructor(factory: PolynomialFactory, coeffs: FiniteField[]) {
    this.factory = factory;
    this.coeffFactory = factory.coeffFactory;

    const normalized = [...coeffs];
    while (normalized.length > 0 && normalized[normalized.length - 1].isZero()) {
      normalized.pop();
    }
    if (normalized.length === 0) {
      this.coeffs = [this.coeffFactory.zero()];
    } else {
      this.coeffs = normalized;
    }
  }

  degree(): number {
    if (this.isZero()) return -1;
    return this.coeffs.length - 1;
  }

  isZero(): boolean {
    return this.coeffs.length === 1 && this.coeffs[0].isZero();
  }

  eq(other: Polynomial): boolean {
    if (this.coeffs.length !== other.coeffs.length) return false;
    for (let i = 0; i < this.coeffs.length; i++) {
      if (!this.coeffs[i].eq(other.coeffs[i])) return false;
    }
    return true;
  }

  clone(): Polynomial {
    return new Polynomial(
      this.factory,
      this.coeffs.map((c) => c.clone())
    );
  }

  add(other: Polynomial): Polynomial {
    const maxLength = Math.max(this.coeffs.length, other.coeffs.length);
    const result: FiniteField[] = [];

    for (let i = 0; i < maxLength; i++) {
      const a = this.coeffs[i] ?? this.coeffFactory.zero();
      const b = other.coeffs[i] ?? this.coeffFactory.zero();
      result.push(a.add(b));
    }

    return new Polynomial(this.factory, result);
  }

  sub(other: Polynomial): Polynomial {
    const maxLength = Math.max(this.coeffs.length, other.coeffs.length);
    const result: FiniteField[] = [];

    for (let i = 0; i < maxLength; i++) {
      const a = this.coeffs[i] ?? this.coeffFactory.zero();
      const b = other.coeffs[i] ?? this.coeffFactory.zero();
      result.push(a.sub(b));
    }

    return new Polynomial(this.factory, result);
  }

  mul(other: Polynomial): Polynomial {
    if (this.isZero() || other.isZero()) {
      return this.factory.zero();
    }
    const resultDegree = this.degree() + other.degree();
    const result = Array(resultDegree + 1).fill(this.coeffFactory.zero());

    for (let i = 0; i < this.coeffs.length; i++) {
      if (this.coeffs[i].isZero()) continue;
      for (let j = 0; j < other.coeffs.length; j++) {
        result[i + j] = result[i + j].add(this.coeffs[i].mul(other.coeffs[j]));
      }
    }

    return new Polynomial(this.factory, result);
  }

  scale(n: bigint): Polynomial {
    return new Polynomial(
      this.factory,
      this.coeffs.map((c) => c.scale(n))
    );
  }

  div(other: Polynomial): Polynomial {
    if (other.isZero()) throw new Error("Division by zero");
    if (this.isZero()) return this.factory.zero();
    if (this.degree() < other.degree()) return this.factory.zero();

    if (other.degree() === 0) {
      const inv = other.coeffs[0].inverse();
      const coeffes = this.coeffs.map((c) => c.mul(inv));
      return new Polynomial(this.factory, coeffes);
    }

    let quotient = this.factory.zero();
    let remainder = this.clone();
    const divisorLeadCoeffInv = other.coeffs[other.degree()].inverse();

    while (!remainder.isZero() && remainder.degree() >= other.degree()) {
      const degreeDiff = remainder.degree() - other.degree();
      const coeff = remainder.coeffs[remainder.degree()].mul(divisorLeadCoeffInv);

      const monomialCoefficients = Array(degreeDiff + 1).fill(this.coeffFactory.zero());
      monomialCoefficients[degreeDiff] = coeff;
      const monomial = new Polynomial(this.factory, monomialCoefficients);

      quotient = quotient.add(monomial);
      remainder = remainder.sub(monomial.mul(other));
    }

    return quotient;
  }

  mod(other: Polynomial): Polynomial {
    if (other.isZero()) throw new Error("Modulo by zero polynomial");
    if (this.isZero()) return this.factory.zero();
    if (this.degree() < other.degree()) return this.clone();

    if (other.degree() === 0) return this.factory.zero();

    let remainder = this.clone();
    const divisorLeadCoeffInv = other.coeffs[other.degree()].inverse();

    while (!remainder.isZero() && remainder.degree() >= other.degree()) {
      const degreeDiff = remainder.degree() - other.degree();
      const coeff = remainder.coeffs[remainder.degree()].mul(divisorLeadCoeffInv);
      const monomialCoefficients = Array(degreeDiff + 1).fill(this.coeffFactory.zero());
      monomialCoefficients[degreeDiff] = coeff;
      const monomial = new Polynomial(this.factory, monomialCoefficients);
      remainder = remainder.sub(monomial.mul(other));
    }

    // ループ終了時の remainder が剰余
    return remainder;
  }

  pow(n: bigint): Polynomial {
    return pow<Polynomial>(this, n);
  }

  eval(x: FiniteField): FiniteField {
    let result = this.coeffFactory.zero();
    for (let i = this.coeffs.length - 1; i >= 0; i--) {
      result = result.mul(x).add(this.coeffs[i]);
    }
    return result;
  }
}
