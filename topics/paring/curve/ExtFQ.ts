import { Polynomial, PolynomialLike } from "./Polynomial";
import { FQ } from "./FQ";

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

  add(a: PolynomialLike, b: PolynomialLike): Polynomial {
    const a_ = Polynomial.mustBePolynomial(a, this.p);
    const b_ = Polynomial.mustBePolynomial(b, this.p);
    return a_.add(b_);
  }
}
