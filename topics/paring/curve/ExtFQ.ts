import { Polynomial, PolynomialLike } from "./Polynomial";
import { FQ } from "./FQ";

export class ExtFQ {
  public readonly baseField: FQ;
  public readonly degree: number;
  public readonly modPoly: Polynomial;

  constructor(baseField: FQ, modPoly: PolynomialLike) {
    const modPoly_ = Polynomial.mustBePolynomial(modPoly, baseField.p);
    if (baseField.p !== modPoly_.p) throw new Error("Must be same field");

    this.baseField = baseField;
    this.degree = modPoly_.degree();
    this.modPoly = modPoly_;
  }

  add(a: PolynomialLike, b: PolynomialLike): Polynomial {
    const a_ = Polynomial.mustBePolynomial(a, this.baseField.p);
    const b_ = Polynomial.mustBePolynomial(b, this.baseField.p);
    return a_.add(b_);
  }
}
