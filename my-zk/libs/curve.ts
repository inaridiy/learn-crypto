import { Field } from "./interface";

export type Point<F extends Field<F, unknown>> = { x: F; y: F };

export class EllipticCurve<F extends Field<F, unknown>> {
  constructor(public readonly a: F, public readonly b: F) {}

  isOnCurve(point: Point<F>): boolean {
    if (point.x.isZero() && point.y.isZero()) return true;
    const left = point.y.pow(2n);
    const right = point.x.pow(3n).add(this.a.mul(point.x)).add(this.b);
    return left.eq(right);
  }

  zero(): EllipticCurvePoint<F> {
    return new EllipticCurvePoint(this, {
      x: this.a.factory.zero(),
      y: this.a.factory.zero(),
    });
  }

  from(value: Point<F>): EllipticCurvePoint<F> {
    return new EllipticCurvePoint(this, value);
  }
}

export class EllipticCurvePoint<F extends Field<F, unknown>> {
  readonly curve: EllipticCurve<F>;
  readonly point: Point<F>;
  readonly f3nf2n: [F, F];

  constructor(curve: EllipticCurve<F>, point: Point<F>) {
    if (!curve.isOnCurve(point)) throw new Error("Point is not on curve");
    this.curve = curve;
    this.point = point;
    this.f3nf2n = [curve.a.factory.one().scale(3n), curve.a.factory.one().scale(2n)];
  }

  isZero(): boolean {
    return this.point.x.isZero() && this.point.y.isZero();
  }

  eq(other: EllipticCurvePoint<F>): boolean {
    return this.point.x.eq(other.point.x) && this.point.y.eq(other.point.y);
  }

  clone(): EllipticCurvePoint<F> {
    return new EllipticCurvePoint(this.curve, this.point);
  }

  add(other: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    if (this.isZero()) return other;
    if (other.isZero()) return this;

    const { x: x1, y: y1 } = this.point;
    const { x: x2, y: y2 } = other.point;
    let x3: F, y3: F;

    if (!x1.eq(x2)) {
      const m = y1.sub(y2).div(x1.sub(x2));
      x3 = m.pow(2n).sub(x1).sub(x2);
      y3 = m.mul(x1.sub(x3)).sub(y1);
    } else if (!y1.eq(y2)) {
      return this.curve.zero();
    } else {
      const [f3n, f2n] = this.f3nf2n; // [3n, 2n];
      const m = x1.pow(2n).mul(f3n).add(this.curve.a).div(y1.mul(f2n));
      x3 = m.pow(2n).sub(x1.mul(f2n));
      y3 = m.mul(x1.sub(x3)).sub(y1);
    }

    return new EllipticCurvePoint(this.curve, { x: x3, y: y3 });
  }

  sub(other: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    if (other.isZero()) return this;
    return this.add(other.neg());
  }

  mul(n: bigint): EllipticCurvePoint<F> {
    if (n < 0n) throw new Error("n must be positive");
    if (n === 0n) return this.curve.zero();

    let result = this.curve.zero();
    let point = this.clone();
    while (n > 0n) {
      if (n % 2n === 1n) result = point.add(result);
      point = point.add(point);
      n >>= 1n;
    }
    return result;
  }

  neg(): EllipticCurvePoint<F> {
    return new EllipticCurvePoint(this.curve, {
      x: this.point.x,
      y: this.point.y.mul(this.curve.a.factory.one().scale(-1n)),
    });
  }
}
