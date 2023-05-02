import { Polynomial, PolynomialLike } from "./Polynomial";
import { extGCD } from "./utils/extendedGCD";
import { FieldFactory, Field } from "./types";
import { fastPow } from "./utils/fastPow";

export type ExtFQLike = ExtFQ | PolynomialLike;

export class ExtFQFactory implements FieldFactory<ExtFQ, ExtFQLike> {
  constructor(public readonly p: bigint, public readonly modPoly: PolynomialLike) {}

  zero(): ExtFQ {
    return ExtFQ.zero(this.p, this.modPoly);
  }

  one(): ExtFQ {
    return new ExtFQ(this.p, this.modPoly, new Polynomial([1n], this.p));
  }

  from(value: PolynomialLike): ExtFQ {
    return new ExtFQ(this.p, this.modPoly, value);
  }
}

export class ExtFQ implements Field<ExtFQ, ExtFQLike> {
  public readonly p: bigint;
  public readonly degree: number;
  public readonly modPoly: Polynomial;

  public readonly value: Polynomial;

  constructor(p: bigint, modPoly: PolynomialLike, value: PolynomialLike) {
    const modPoly_ = Polynomial.mustBePolynomial(modPoly, p);
    const value_ = Polynomial.mustBePolynomial(value, p);
    if (p !== modPoly_.p) throw new Error("Must be same field");
    if (p !== value_.p) throw new Error("Must be same field");

    this.p = p;
    this.degree = modPoly_.degree();
    this.modPoly = modPoly_;
    this.value = value_.mod(modPoly_);
  }

  static mustBeExtFQ(other: ExtFQLike, p: bigint, modPoly: PolynomialLike): ExtFQ {
    if (other instanceof ExtFQ) return other;
    else return new ExtFQ(p, modPoly, other);
  }

  static zero = (p: bigint, modPoly: PolynomialLike): ExtFQ =>
    new ExtFQ(p, modPoly, Polynomial.zero(p));
  static one = (p: bigint, modPoly: PolynomialLike): ExtFQ =>
    new ExtFQ(p, modPoly, Polynomial.one(p));
  zero = (): ExtFQ => ExtFQ.zero(this.p, this.modPoly);
  one = (): ExtFQ => ExtFQ.one(this.p, this.modPoly);

  extend(value: PolynomialLike): ExtFQ {
    const value_ = Polynomial.mustBePolynomial(value, this.p);
    if (this.p !== value_.p) throw new Error("Must be same field");
    return new ExtFQ(this.p, this.modPoly, value_);
  }

  clone(): ExtFQ {
    return new ExtFQ(this.p, this.modPoly, this.value);
  }

  eq(other: ExtFQLike): boolean {
    const other_ = ExtFQ.mustBeExtFQ(other, this.p, this.modPoly);
    return this.value.eq(other_.value);
  }

  isZero(): boolean {
    return this.value.isZero();
  }

  /**
   *
   * @returns
   */
  inverse(): ExtFQ {
    const [g, x, _] = extGCD(this.value, this.modPoly);
    //公約数が定数項の場合は逆元が存在しない
    if (g.degree() !== 0) throw new Error("Not invertible");
    //定数項が1でない場合は、定数項を1にするために定数項で割って正規化する?
    else if (g.coefficients[0].n !== 1n) return this.extend(x.div(g));
    //定数項が1の場合はそのまま返す
    else return this.extend(x);
  }

  mod(): ExtFQ {
    return this.extend(this.value.mod(this.modPoly));
  }

  add(other: ExtFQLike): ExtFQ {
    const other_ = ExtFQ.mustBeExtFQ(other, this.p, this.modPoly);
    const result = this.value.add(other_.value);
    return this.extend(result);
  }

  sub(other: ExtFQLike): ExtFQ {
    const other_ = ExtFQ.mustBeExtFQ(other, this.p, this.modPoly);
    const result = this.value.sub(other_.value);
    return this.extend(result);
  }

  mul(other: ExtFQLike): ExtFQ {
    const other_ = ExtFQ.mustBeExtFQ(other, this.p, this.modPoly);
    const result = this.value.mul(other_.value);
    return this.extend(result);
  }

  div(other: ExtFQLike): ExtFQ {
    const other_ = ExtFQ.mustBeExtFQ(other, this.p, this.modPoly);
    const result = this.value.mul(other_.inverse().value);
    return this.extend(result);
  }

  pow(n: bigint): ExtFQ {
    return fastPow(this, n);
  }

  toString(): string {
    return this.value.toString();
  }
}
