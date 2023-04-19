import { extendedGCD } from "./extendedGCD";

export class Finite {
  private _inverseCache: Map<bigint, bigint> = new Map();

  constructor(readonly prime: bigint) {}

  public inverse(n: bigint): bigint {
    if (this._inverseCache.has(n)) return this._inverseCache.get(n)!;
    const [gcd, x] = extendedGCD(n, this.prime);
    if (gcd !== 1n) throw new Error("No inverse");
    const result = this.mod(x);
    this._inverseCache.set(x, result);
    return result;
  }

  public mod(x: bigint): bigint {
    const r = x % this.prime;
    return r < 0n ? r + this.prime : r;
  }

  public add(x: bigint, y: bigint): bigint {
    return this.mod(x + y);
  }

  public mul(x: bigint, y: bigint): bigint {
    return this.mod(x * y);
  }

  public div(x: bigint, y: bigint): bigint {
    return this.mod(x * this.inverse(y));
  }
}
