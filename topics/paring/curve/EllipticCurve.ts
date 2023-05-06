import { ExtFQ, ExtFQLike } from "./ExtFQ";
import { FieldFactory } from "./types";

export type PointCoord = Readonly<{
  x: ExtFQ;
  y: ExtFQ;
}>;

export class EllipticCurve {
  readonly a: ExtFQ;
  readonly b: ExtFQ;
  readonly fq: FieldFactory<ExtFQ, ExtFQLike>;

  zero: () => PointCoord = () => ({ x: this.fq.zero(), y: this.fq.zero() });

  constructor(a: ExtFQLike, b: ExtFQLike, fieldFactory: FieldFactory<ExtFQ, ExtFQLike>) {
    this.a = fieldFactory.from(a);
    this.b = fieldFactory.from(b);
    this.fq = fieldFactory;
  }

  isOnCurve(point: PointCoord): boolean {
    if (point.x.isZero() && point.y.isZero()) return true;
    const left = point.y.pow(2n);
    const right = point.x.pow(3n).add(this.a.mul(point.x)).add(this.b);
    return left.eq(right);
  }

  add(p1: PointCoord, p2: PointCoord): PointCoord {
    if (p1.x.isZero() && p1.y.isZero()) return p2;
    if (p2.x.isZero() && p2.y.isZero()) return p1;
    if (!this.isOnCurve(p1) || !this.isOnCurve(p2)) throw new Error("Points must be on curve");

    const { x: x1, y: y1 } = p1;
    const { x: x2, y: y2 } = p2;
    let x3: ExtFQ, y3: ExtFQ;

    if (!x1.eq(x2)) {
      const m = y1.sub(y2).div(x1.sub(x2));
      x3 = m.pow(2n).sub(x1).sub(x2);
      y3 = m.mul(x1.sub(x3)).sub(y1);
    } else if (!y1.eq(y2)) {
      return this.zero();
    } else {
      const [f3n, f2n] = [this.fq.from(3n), this.fq.from(2n)]; // [3n, 2n];
      const m = x1.pow(2n).mul(f3n).add(this.a).div(y1.mul(f2n));
      x3 = m.pow(2n).sub(x1.mul(f2n));
      y3 = m.mul(x1.sub(x3)).sub(y1);
    }

    return { x: x3, y: y3 };
  }

  sub(p1: PointCoord, p2: PointCoord): PointCoord {
    if (p2.x.isZero() && p2.x.isZero()) return p1;
    return this.add(p1, {
      x: p2.x,
      y: p2.y.mul(this.fq.from(-1n)), // -p2.y
    });
  }

  mul(point: PointCoord, n: bigint) {
    if (n < 0n) throw new Error("n must be positive");
    if (n === 0n) return this.zero();

    let result = this.zero();
    while (n > 0n) {
      if (n % 2n === 1n) result = this.add(result, point);
      point = this.add(point, point);
      n >>= 1n;
    }
    return result;
  }
}

export class CurvePoint {
  x: ExtFQ;
  y: ExtFQ;
  curve: EllipticCurve;

  constructor(x: ExtFQLike, y: ExtFQLike, curve: EllipticCurve) {
    this.x = curve.fq.from(x);
    this.y = curve.fq.from(y);
    this.curve = curve;
  }

  isOnCurve(): boolean {
    return this.curve.isOnCurve(this);
  }

  add(p: CurvePoint): CurvePoint {
    const { x, y } = this.curve.add(this, p);
    return new CurvePoint(x, y, this.curve);
  }

  sub(p: CurvePoint): CurvePoint {
    const { x, y } = this.curve.sub(this, p);
    return new CurvePoint(x, y, this.curve);
  }

  mul(n: bigint): CurvePoint {
    const { x, y } = this.curve.mul(this, n);
    return new CurvePoint(x, y, this.curve);
  }
}
