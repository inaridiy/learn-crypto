import { extendedGCD } from "./extendedGCD";

type Point = Readonly<{
  x: bigint;
  y: bigint;
}>;

export class EllipticCurve {
  private a: bigint;
  private b: bigint;
  private p: bigint;
  private _inverseCache: Map<bigint, bigint> = new Map();
  private _onCurveCache: Map<bigint, boolean> = new Map();

  static BIG_PRIME = 2n ** 256n - 2n ** 32n - 977n;
  static ZERO_POINT: Point = { x: 0n, y: 0n };

  static SECP256K1: EllipticCurve = new EllipticCurve(0n, 7n);
  static SECP256K1_G: Point = {
    x: 55066263022277343669578718895168534326250603453777594175500187360389116729240n,
    y: 32670510020758816978083085130507043184471273380659243275938904335757337482424n,
  };

  constructor(a: bigint, b: bigint, p: bigint = EllipticCurve.BIG_PRIME) {
    this.a = a;
    this.b = b;
    this.p = p;
  }

  private mod(x: bigint): bigint {
    const r = x % this.p;
    return r < 0n ? r + this.p : r;
  }

  private inverse(n: bigint): bigint {
    if (this._inverseCache.has(n)) return this._inverseCache.get(n)!;
    const [gcd, x] = extendedGCD(n, this.p);
    if (gcd !== 1n) throw new Error("No inverse");
    const result = this.mod(x);
    this._inverseCache.set(x, result);
    return result;
  }

  public isOnCurve(point: Point): boolean {
    if (this._onCurveCache.has(point.x)) return this._onCurveCache.get(point.x)!;
    const left = this.mod(point.y ** 2n);
    const right = this.mod(point.x ** 3n + this.a * point.x + this.b);
    const result = left === right;
    this._onCurveCache.set(point.x, result);
    return result;
  }

  private _add(p1: Point, p2: Point): Point {
    if (p1.x === 0n && p1.y === 0n) return p2;
    if (p2.x === 0n && p2.y === 0n) return p1;
    if (!this.isOnCurve(p1) || !this.isOnCurve(p2)) throw new Error("Points must be on curve");
    if (p1.x > p2.x) [p1, p2] = [p2, p1]; //なんかこうしないと通らない

    let x3: bigint, y3: bigint;

    if (p1.x !== p2.x) {
      const m = this.mod((p2.y - p1.y) * this.inverse(p2.x - p1.x));
      x3 = this.mod(m ** 2n - p1.x - p2.x);
      y3 = this.mod(m * (p1.x - x3) - p1.y);
    } else {
      const m = this.mod((3n * p1.x * p1.x + this.a) * this.inverse(2n * p1.y));
      x3 = this.mod(m * m - 2n * p1.x);
      y3 = this.mod(m * (p1.x - x3) - p1.y);
    }

    return { x: x3, y: y3 };
  }

  public add(...points: Point[]): Point {
    if (points.length === 0) return EllipticCurve.ZERO_POINT;
    if (points.length === 1) return points[0];
    let point = points[0];
    for (let i = 1; i < points.length; i++) {
      point = this._add(point, points[i]);
    }
    return point;
  }

  public multiply(point: Point, n: bigint): Point {
    if (n < 0n) throw new Error("n must be positive");
    if (n === 0n) return EllipticCurve.ZERO_POINT;

    let result = EllipticCurve.ZERO_POINT;
    while (n > 0n) {
      if (n % 2n === 1n) result = this._add(result, point);
      point = this._add(point, point);
      n >>= 1n;
    }
    return result;
  }
}
