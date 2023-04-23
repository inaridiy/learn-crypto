import { extendedGCD } from "./extendedGCD";

export class FiniteField {
  constructor(readonly prime: bigint) {}

  public inverse(n: bigint): bigint {
    const [gcd, x] = extendedGCD(n, this.prime);
    if (gcd !== 1n) throw new Error("No inverse");
    const result = this.mod(x);
    return result;
  }

  public mod(x: bigint): bigint {
    const r = x % this.prime;
    return r < 0n ? r + this.prime : r;
  }

  public add(x: bigint, y: bigint): bigint {
    return this.mod(x + y);
  }

  public sub(x: bigint, y: bigint): bigint {
    return this.mod(x - y);
  }

  public mul(x: bigint, y: bigint): bigint {
    return this.mod(x * y);
  }

  public div(x: bigint, y: bigint): bigint {
    return this.mod(x * this.inverse(y));
  }
}
