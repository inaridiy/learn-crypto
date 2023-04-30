import { Field } from "./Field";

export type Point = Readonly<{
  x: bigint;
  y: bigint;
}>;

export class EllipticCurve {
  readonly a: bigint;
  readonly b: bigint;
  readonly p: bigint;
  readonly fp: Field;
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
    this.fp = new Field(p);
  }

  public isOnCurve(point: Point): boolean {
    if (point.x === 0n && point.y === 0n) return true;
    if (this._onCurveCache.has(point.x)) return this._onCurveCache.get(point.x)!;
    const left = this.fp.mod(point.y ** 2n);
    const right = this.fp.mod(point.x ** 3n + this.a * point.x + this.b);
    const result = left === right;
    this._onCurveCache.set(point.x, result);
    return result;
  }

  private _add(p1: Point, p2: Point): Point {
    if (p1.x === 0n && p1.y === 0n) return p2;
    if (p2.x === 0n && p2.y === 0n) return p1;
    if (!this.isOnCurve(p1) || !this.isOnCurve(p2)) throw new Error("Points must be on curve");
    if (p1.x > p2.x) [p1, p2] = [p2, p1]; //巡回群内での演算の高速化をしつつ、アンダーフローを防止する。

    let x3: bigint, y3: bigint;

    if (p1.x !== p2.x) {
      const m = this.fp.div(p2.y - p1.y, p2.x - p1.x);
      x3 = this.fp.mod(m ** 2n - p1.x - p2.x);
      y3 = this.fp.mod(m * (p1.x - x3) - p1.y);
    } else if (p1.y !== p2.y) {
      return EllipticCurve.ZERO_POINT;
    } else {
      const m = this.fp.div(3n * p1.x * p1.x + this.a, 2n * p1.y);
      x3 = this.fp.mod(m * m - 2n * p1.x);
      y3 = this.fp.mod(m * (p1.x - x3) - p1.y);
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

  public sub(p1: Point, p2: Point): Point {
    if (p2.x === 0n && p2.y === 0n) return p1;
    return this._add(p1, { x: p2.x, y: this.fp.mod(-p2.y) });
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
