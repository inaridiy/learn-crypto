import { Polynomial, PolynomialLike } from "./Polynomial";
import { extGCD } from "./extendedGCD";

export class ExtFQ {
  public readonly p: bigint;
  public readonly degree: number;
  public readonly modPoly: Polynomial;

  constructor(p: bigint, modPoly: PolynomialLike) {
    const modPoly_ = Polynomial.mustBePolynomial(modPoly, p);
    if (p !== modPoly_.p) throw new Error("Must be same field");

    this.p = p;
    this.degree = modPoly_.degree();
    this.modPoly = modPoly_;
  }

  inverse(a: PolynomialLike): Polynomial {
    const a_ = Polynomial.mustBePolynomial(a, this.p);
    const [gcd, x] = extGCD(a_, this.modPoly);
    if (gcd.degree() !== 0) throw new Error("No inverse");
    return x.mod(this.modPoly);
  }

  mod(a: PolynomialLike): Polynomial {
    const a_ = Polynomial.mustBePolynomial(a, this.p);
    return a_.mod(this.modPoly);
  }

  private _normalize(a: PolynomialLike, b: PolynomialLike): [Polynomial, Polynomial] {
    const a_ = Polynomial.mustBePolynomial(a, this.p);
    const b_ = Polynomial.mustBePolynomial(b, this.p);
    return [a_, b_];
  }

  add(a: PolynomialLike, b: PolynomialLike): Polynomial {
    const [a_, b_] = this._normalize(a, b);
    return a_.add(b_);
  }

  sub(a: PolynomialLike, b: PolynomialLike): Polynomial {
    const [a_, b_] = this._normalize(a, b);
    return a_.sub(b_);
  }

  mul(a: PolynomialLike, b: PolynomialLike): Polynomial {
    const [a_, b_] = this._normalize(a, b);
    return a_.mul(b_).mod(this.modPoly);
  }

  div(a: PolynomialLike, b: PolynomialLike): Polynomial {
    const [a_, b_] = this._normalize(a, b);
    return a_.mul(this.inverse(b_)).mod(this.modPoly);
  }
}
