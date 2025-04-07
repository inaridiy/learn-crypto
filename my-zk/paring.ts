import { BLS12_381_FQ12, CURVE_ORDER, FIELD_MODULUS, FQ12, POLY } from "./bls12-381";
import { EllipticCurvePoint } from "./curve";
import { ExtendedFiniteField } from "./extended-finite-field";

const ATE_LOOP_COUNT = 15132376222941642752n;
const LOG_ATE_LOOP_COUNT = 62n;

const embedFQ12 = (P: EllipticCurvePoint<ExtendedFiniteField>) => {
  return BLS12_381_FQ12.from({
    x: FQ12.from(P.point.x.n),
    y: FQ12.from(P.point.y.n),
  });
};

const twist = (P: EllipticCurvePoint<ExtendedFiniteField>) => {
  const xc = [P.point.x.n.coeffs[0].sub(P.point.x.n.coeffs[1]), P.point.x.n.coeffs[1]];
  const yc = [P.point.y.n.coeffs[0].sub(P.point.y.n.coeffs[1]), P.point.y.n.coeffs[1]];

  const nx = FQ12.from(POLY.from([xc[0].n, 0n, 0n, 0n, 0n, 0n, xc[1].n, 0n, 0n, 0n, 0n, 0n, 0n]));
  const ny = FQ12.from(POLY.from([yc[0].n, 0n, 0n, 0n, 0n, 0n, yc[1].n, 0n, 0n, 0n, 0n, 0n, 0n]));

  const mx = nx.div(FQ12.from(POLY.from([0n, 0n, 1n])));
  const my = ny.div(FQ12.from(POLY.from([0n, 0n, 0n, 1n])));

  return BLS12_381_FQ12.from({ x: mx, y: my });
};

const lineFunction = (
  P1: EllipticCurvePoint<ExtendedFiniteField>,
  P2: EllipticCurvePoint<ExtendedFiniteField>,
  Q: EllipticCurvePoint<ExtendedFiniteField>
) => {
  if (P1.point.x.n.isZero() || P2.point.x.n.isZero() || Q.point.x.n.isZero()) return FQ12.one();

  const { x: x1, y: y1 } = P1.point;
  const { x: x2, y: y2 } = P2.point;
  const { x: xq, y: yq } = Q.point;

  if (!x1.eq(x2)) {
    const m = y2.sub(y1).div(x2.sub(x1));
    return xq.sub(x1).mul(m).add(y1).sub(yq);
  } else if (y1.eq(y2)) {
    const m = x1.pow(2n).scale(3n).div(y1.scale(2n));
    return xq.sub(x1).mul(m).add(y1).sub(yq);
  } else {
    return xq.sub(x1);
  }
};

export const millerLoop = (
  P: EllipticCurvePoint<ExtendedFiniteField>,
  Q: EllipticCurvePoint<ExtendedFiniteField>
) => {
  let R = Q;
  let f = FQ12.one();

  for (let i = LOG_ATE_LOOP_COUNT; i >= 0n; i--) {
    f = f.pow(2n).mul(lineFunction(R, R, P));
    R = R.mul(2n);

    if (ATE_LOOP_COUNT & (1n << i)) {
      f = f.mul(lineFunction(R, Q, P));
      R = R.add(Q);
    }
  }

  return f.pow((FIELD_MODULUS ** 12n - 1n) / CURVE_ORDER);
};

export const pairing = (
  P: EllipticCurvePoint<ExtendedFiniteField>,
  Q: EllipticCurvePoint<ExtendedFiniteField>
) => {
  const P_ = P.curve.a.mod.degree() === 2 ? twist(P) : embedFQ12(P);
  const Q_ = Q.curve.a.mod.degree() === 2 ? twist(Q) : embedFQ12(Q);

  return millerLoop(P_, Q_);
};
