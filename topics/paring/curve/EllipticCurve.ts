import { Field, FieldFactory } from "./types";

export type PointCoord<T extends Field<any, any> = Field<any, any>> = Readonly<{
  x: T;
  y: T;
}>;

export class EllipticCurve<T, TLike extends T> {
  readonly a: Field<T, TLike>;
  readonly b: Field<T, TLike>;
  readonly fq: FieldFactory<T, TLike>;

  zero: () => PointCoord<Field<T, TLike>> = () => ({ x: this.fq.zero(), y: this.fq.zero() });

  constructor(a: TLike, b: TLike, fieldFactory: FieldFactory<T, TLike>) {
    this.a = fieldFactory.from(a);
    this.b = fieldFactory.from(b);
    this.fq = fieldFactory;
  }

  isOnCurve(point: PointCoord<Field<T, TLike>>): boolean {
    if (point.x.isZero() && point.y.isZero()) return true;
    const left = point.y.pow(2n);
    const right = point.x.pow(3n).add(this.a.mul(point.x)).add(this.b);
    return left.eq(right);
  }

  add(
    p1: PointCoord<Field<T, TLike>>,
    p2: PointCoord<Field<T, TLike>>
  ): PointCoord<Field<T, TLike>> {
    if (p1.x.isZero() && p1.y.isZero()) return p2;
    if (p2.x.isZero() && p2.y.isZero()) return p1;
    if (!this.isOnCurve(p1) || !this.isOnCurve(p2)) throw new Error("Points must be on curve");

    const { x: x1, y: y1 } = p1;
    const { x: x2, y: y2 } = p2;
    let x3: Field<T, TLike>, y3: Field<T, TLike>;

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
}
