import { extGCD } from "./utils/extendedGCD";
import { fastPow } from "./utils/fastPow";

export type FQLike = FQ | bigint;

export class FQFactory {
  constructor(public readonly p: bigint) { }

  zero(): FQ {
    return FQ.zero(this.p);
  }

  one(): FQ {
    return new FQ(this.p, 1n);
  }

  from(value: FQLike): FQ {
    return FQ.mustBeFQ(value, this.p);
  }
}

export class FQ {
  public readonly p: bigint;
  public readonly n: bigint;

  constructor(p: bigint, n: bigint) {
    this.p = p;
    this.n = this._mod(n);
  }

  static mustBeFQ(other: FQ | bigint, p: bigint): FQ {
    if (other instanceof FQ)
      if (p !== other.p) throw new Error("FQs must have same p");
      else return other;
    else return new FQ(p, other);
  }
  mustBeFQ(other: FQ | bigint): FQ {
    return FQ.mustBeFQ(other, this.p);
  }

  static zero = (p: bigint): FQ => new FQ(p, 0n);
  static one = (p: bigint): FQ => new FQ(p, 1n);
  zero = (): FQ => FQ.zero(this.p);
  one = (): FQ => FQ.one(this.p);

  extend(other: FQ | bigint): FQ {
    const other_ = this.mustBeFQ(other);
    return new FQ(this.p, other_.n);
  }

  isZero(): boolean {
    return this.n === 0n;
  }

  clone(): FQ {
    return new FQ(this.p, this.n);
  }

  eq(other: FQ | bigint): boolean {
    const other_ = this.mustBeFQ(other);
    return this.n === other_.n;
  }

  private _mod(n: bigint): bigint {
    const r = n % this.p;
    return 0n > r ? r + this.p : r;
  }

  mod(): FQ {
    return new FQ(this.p, this._mod(this.n));
  }

  private _inverse(): bigint {
    const [gcd, x] = extGCD(this.n, this.p);
    if (gcd !== 1n) throw new Error("No inverse");
    const result = this._mod(x);
    return result;
  }

  inverse(): FQ {
    return new FQ(this.p, this._inverse());
  }

  add(other: FQ | bigint): FQ {
    const other_ = this.mustBeFQ(other);
    return new FQ(this.p, this._mod(this.n + other_.n));
  }

  sub(other: FQ | bigint): FQ {
    const other_ = this.mustBeFQ(other);
    return new FQ(this.p, this._mod(this.n - other_.n));
  }

  mul(other: FQ | bigint): FQ {
    const other_ = this.mustBeFQ(other);
    return new FQ(this.p, this._mod(this.n * other_.n));
  }

  div(other: FQ | bigint): FQ {
    const other_ = this.mustBeFQ(other);
    return new FQ(this.p, this._mod(this.n * other_.inverse().n));
  }

  pow(n: bigint): FQ {

    return fastPow(this, n);
  }
}
