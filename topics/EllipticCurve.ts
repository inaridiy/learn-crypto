type Point = Readonly<{
  x: bigint;
  y: bigint;
}>;

export class EllipticCurve {
  private a: bigint;
  private b: bigint;
  private p: bigint;
  private _modCache: Map<bigint, bigint> = new Map();
  private _inverseCache: Map<bigint, bigint> = new Map();

  static BIG_PRIME = 2n ** 256n - 2n ** 32n - 977n;
  static ZERO_POINT: Point = { x: 0n, y: 0n };

  constructor(a: bigint, b: bigint, p: bigint = EllipticCurve.BIG_PRIME) {
    this.a = a;
    this.b = b;
    this.p = p;
  }

  //負の値の剰余を正の値にする
  private mod(n: bigint): bigint {
    if (this._modCache.has(n)) return this._modCache.get(n)!;
    else return ((n % this.p) + this.p) % this.p;
  }

  private inverse(n: bigint): bigint {
    if (this._inverseCache.has(n)) return this._inverseCache.get(n)!;
    let [g, x, _] = this.extendedGCD(n, this.p);
    if (g !== 1n) throw new Error("Inverse does not exist.");

    const result = this.mod(x);
    this._inverseCache.set(n, result);
    return result;
  }

  private extendedGCD(a: bigint, b: bigint): [bigint, bigint, bigint] {
    let x = 0n;
    let y = 1n;
    let u = 1n;
    let v = 0n;

    while (a !== 0n) {
      const q = b / a;
      const r = b % a;
      const m = x - u * q;
      const n = y - v * q;

      b = a;
      a = r;
      x = u;
      y = v;
      u = m;
      v = n;
    }

    return [b, x, y];
  }

  public isOnCurve(point: Point): boolean {
    const { x, y } = point;
    const left = this.mod(y * y);
    const right = this.mod(x * x * x + this.a * x + this.b);
    return left === right;
  }

  private _add(point1: Point, point2: Point): Point {
    if (!this.isOnCurve(point1) || !this.isOnCurve(point2)) {
      throw new Error("The points are not on the curve.");
    }

    let x3, y3: bigint;

    if (point1.x === point2.x && point1.y === point2.y) {
      if (point1.y === 0n) return EllipticCurve.ZERO_POINT;

      const m = this.mod((3n * point1.x * point1.x + this.a) * this.inverse(2n * point1.y));
      x3 = this.mod(m * m - 2n * point1.x);
      y3 = this.mod(m * (point1.x - x3) - point1.y);
    } else {
      if (point1.x === point2.x) return EllipticCurve.ZERO_POINT;

      const m = this.mod((point2.y - point1.y) * this.inverse(point2.x - point1.x));
      x3 = this.mod(m * m - point1.x - point2.x);
      y3 = this.mod(m * (point1.x - x3) - point1.y);
    }

    return { x: x3, y: y3 };
  }

  public add(...points: Point[]): Point {
    if (points.length === 0) return EllipticCurve.ZERO_POINT;
    let result = points[0];
    for (const point of points.slice(1)) result = this._add(result, point);
    return result;
  }

  public multiply(point: Point, k: bigint): Point {
    if (!this.isOnCurve(point)) throw new Error("The point is not on the curve.");
    if (k === 0n) return EllipticCurve.ZERO_POINT;

    let result: Point = { x: point.x, y: point.y };
    let tempPoint: Point = { x: point.x, y: point.y };

    k -= 1n;

    while (k > 0n) {
      if (k % 2n === 1n) {
        result = this._add(result, tempPoint);
      }
      tempPoint = this._add(tempPoint, tempPoint);
      k /= 2n;
    }

    return result || EllipticCurve.ZERO_POINT;
  }
}
