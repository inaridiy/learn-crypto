import { extendedGCD } from "./extendedGCD";

export class FQ {
  public readonly p: bigint;
  public readonly n: bigint;

  static zero(p: bigint): FQ {
    return new FQ(p, 0n);
  }

  constructor(p: bigint, n: bigint) {
    this.p = p;
    this.n = n;
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

  isZero(): boolean {
    return this.n === 0n;
  }

  clone(): FQ {
    return new FQ(this.p, this.n);
  }

  private _mod(n: bigint): bigint {
    const r = n % this.p;
    return 0n > r ? r + this.p : r;
  }

  mod(other: FQ | bigint): FQ {
    const other_ = this.mustBeFQ(other);
    return new FQ(this.p, this._mod(other_.n));
  }

  private _inverse(): bigint {
    const [gcd, x] = extendedGCD(this.n, this.p);
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
    if (8n > n) return this.mod(this.n ** n); //nが小さい場合は普通に計算する

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
