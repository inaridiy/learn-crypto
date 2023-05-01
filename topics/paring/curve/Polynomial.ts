import { FQ } from "./FQ";

export type PolynomialLike = Polynomial | bigint[] | FQ[];

export class Polynomial {
  public readonly p: bigint;
  public readonly coefficients: FQ[];

  constructor(coefficients_: bigint[] | FQ[], p: bigint) {
    this.p = p;
    this.coefficients = Polynomial.normalize(coefficients_, p);
  }

  static normalize(coefficients: bigint[] | FQ[], p: bigint): FQ[] {
    const coefficients_ = coefficients.map((c) => FQ.mustBeFQ(c, p));
    let result = [...coefficients_];
    while (result.length > 1 && result[result.length - 1].isZero()) result.pop();
    return result;
  }

  static mustBePolynomial(polynomial: PolynomialLike, p: bigint): Polynomial {
    if (polynomial instanceof Polynomial) return polynomial;
    else if (Array.isArray(polynomial)) return new Polynomial(polynomial, p);
    else throw new Error("Not a polynomial");
  }

  mustBePolynomial(other: PolynomialLike): Polynomial {
    return Polynomial.mustBePolynomial(other, this.p);
  }

  clone(): Polynomial {
    return new Polynomial(this.coefficients, this.p);
  }

  degree(): number {
    return this.coefficients.length - 1;
  }

  mod(other: PolynomialLike): Polynomial {
    const other_ = this.mustBePolynomial(other);
    if (this.degree() < other_.degree()) return this;

    let result = this.clone();
    while (result.degree() >= other_.degree()) {
      const degreeDiff = result.degree() - other_.degree();
      const coeff = result.coefficients[result.degree()].div(other_.coefficients[other_.degree()]);
      const monomialCoefficients = Array(degreeDiff + 1).fill(FQ.zero(this.p));
      monomialCoefficients[degreeDiff] = coeff;
      const monomial = new Polynomial(monomialCoefficients, this.p);
      result = result.sub(monomial.mul(other_));
    }

    return result;
  }

  add(other: PolynomialLike): Polynomial {
    const other_ = this.mustBePolynomial(other);
    const maxLen = Math.max(this.coefficients.length, other_.coefficients.length);
    const result: FQ[] = Array(maxLen).fill(FQ.zero(this.p));

    for (let i = 0; i < maxLen; i++) {
      const a = this.coefficients[i] || FQ.zero(this.p);
      const b = other_.coefficients[i] || FQ.zero(this.p);
      result[i] = a.add(b);
    }

    return new Polynomial(result, this.p);
  }

  sub(other: PolynomialLike): Polynomial {
    const other_ = this.mustBePolynomial(other);
    const maxLen = Math.max(this.coefficients.length, other_.coefficients.length);
    const result = Array(maxLen).fill(FQ.zero(this.p));

    for (let i = 0; i < maxLen; i++) {
      const a = this.coefficients[i] || FQ.zero(this.p);
      const b = other_.coefficients[i] || FQ.zero(this.p);
      result[i] = a.sub(b);
    }

    return new Polynomial(result, this.p);
  }

  mul(other: PolynomialLike): Polynomial {
    const other_ = this.mustBePolynomial(other);
    const resultDegree = this.coefficients.length + other_.coefficients.length - 1;
    const result = Array(resultDegree).fill(FQ.zero(this.p));

    for (let i = 0; i < this.coefficients.length; i++) {
      for (let j = 0; j < other_.coefficients.length; j++) {
        result[i + j] = result[i + j].add(this.coefficients[i].mul(other_.coefficients[j]));
      }
    }

    return new Polynomial(result, this.p);
  }

  eval(x: FQ | bigint): FQ {
    const x_ = FQ.mustBeFQ(x, this.p);
    let result = FQ.zero(this.p);
    for (let i = this.coefficients.length - 1; i >= 0; i--) {
      result = result.mul(x_).add(this.coefficients[i]);
    }
    return result;
  }

  toString(): string {
    const terms = this.coefficients.map((c, i) => {
      if (c.isZero())
        if (this.coefficients.length === 1) return "0";
        else return "";
      else if (i === 0) return `${c.n}`;
      else if (i === 1) return `${c.n === 1n ? "" : c.n}x`;
      else return `${c.n === 1n ? "" : c.n}x^${i}`;
    });
    return terms.reverse().join("+").replace(/\+\-/g, "-").replace(/\+$/g, "").replace(/\++/g, "+");
  }
}
