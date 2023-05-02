import { Field, FieldFactory } from "./types";
import { extGCD } from "./extendedGCD";

export type FQLike = FQ | bigint;

export class FQFactory implements FieldFactory<FQ, FQLike> {
  constructor(public readonly p: bigint) {}

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

export class FQ implements Field<FQ, FQLike> {
  public readonly p: bigint;
  public readonly n: bigint;

  static zero(p: bigint): FQ {
    return new FQ(p, 0n);
  }

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

  extend(other: FQ | bigint): FQ {
    const other_ = this.mustBeFQ(other);
    return new FQ(this.p, other_.n);
  }

  mustBeFQ(other: FQ | bigint): FQ {
    return FQ.mustBeFQ(other, this.p);
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
    if (n === 0n) return new FQ(this.p, 1n);
    if (n === 1n) return this.clone();
    if (8n > n) return this.extend(this._mod(this.n ** n)); //nが小さい場合は普通に計算する

    let result = this.clone();
    let m = n;
    while (m > 1n) {
      if (m % 2n === 1n) {
        result = result.mul(this);
        m -= 1n;
      }
      result = result.mul(result);
      m /= 2n;
    }
    return result;
  }
}
