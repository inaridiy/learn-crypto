import {
  BLS12_381_FQ12,
  BLS12_381_G1,
  BLS12_381_G2,
  CURVE_ORDER,
  FIELD_MODULUS,
  FQ12,
  POLY,
} from "./curves/bls12-381";
import { ECPointOnExtendedFF } from "./elliptic-curve";

const ATE_LOOP_COUNT = 15132376222941642752n;
const LOG_ATE_LOOP_COUNT = 62n;

const embedFQ12 = (P: ECPointOnExtendedFF) => {
  return BLS12_381_FQ12.from({
    x: FQ12.from(P.x.n),
    y: FQ12.from(P.y.n),
  });
};

const twist = (P: ECPointOnExtendedFF) => {
  const xZero = P.x.n.structure.coeffField.zero();
  const yZero = P.y.n.structure.coeffField.zero();

  const [x0, x1] = [P.x.n.coeffs[0] ?? xZero, P.x.n.coeffs[1] ?? xZero];
  const [y0, y1] = [P.y.n.coeffs[0] ?? yZero, P.y.n.coeffs[1] ?? yZero];

  const xc = [x0.sub(x1), x1];
  const yc = [y0.sub(y1), y1];

  const nx = FQ12.from(POLY.from([xc[0].n, 0n, 0n, 0n, 0n, 0n, xc[1].n, 0n, 0n, 0n, 0n, 0n, 0n]));
  const ny = FQ12.from(POLY.from([yc[0].n, 0n, 0n, 0n, 0n, 0n, yc[1].n, 0n, 0n, 0n, 0n, 0n, 0n]));

  const mx = nx.div(FQ12.from(POLY.from([0n, 0n, 1n])));
  const my = ny.div(FQ12.from(POLY.from([0n, 0n, 0n, 1n])));

  return BLS12_381_FQ12.from({ x: mx, y: my });
};

const lineFunction = (P1: ECPointOnExtendedFF, P2: ECPointOnExtendedFF, Q: ECPointOnExtendedFF) => {
  if (P1.x.n.isZero() || P2.x.n.isZero() || Q.x.n.isZero()) return FQ12.one();

  const { x: x1, y: y1 } = P1;
  const { x: x2, y: y2 } = P2;
  const { x: xq, y: yq } = Q;

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

export const millerLoop = (P: ECPointOnExtendedFF, Q: ECPointOnExtendedFF) => {
  let R = Q;
  let f = FQ12.one();

  for (let i = LOG_ATE_LOOP_COUNT; i >= 0n; i--) {
    f = f.pow(2n).mul(lineFunction(R, R, P));
    R = R.scale(2n);

    if (ATE_LOOP_COUNT & (1n << i)) {
      f = f.mul(lineFunction(R, Q, P));
      R = R.add(Q);
    }
  }

  return f.pow((FIELD_MODULUS ** 12n - 1n) / CURVE_ORDER);
};

export const pairing = (P: ECPointOnExtendedFF, Q: ECPointOnExtendedFF) => {
  const P_ = P.structure.a.p.degree() === 2 ? twist(P) : embedFQ12(P);
  const Q_ = Q.structure.a.p.degree() === 2 ? twist(Q) : embedFQ12(Q);

  return millerLoop(P_, Q_);
};

if (import.meta.vitest) {
  const { test, expect } = import.meta.vitest;
  test("pairing", () => {
    const a = 123n;
    const b = 456n;

    const aP = BLS12_381_G1.generator().scale(a);
    const bQ = BLS12_381_G2.generator().scale(b);

    const e1 = pairing(BLS12_381_G1.generator(), BLS12_381_G2.generator());
    const e2 = pairing(aP, bQ);

    expect(e1.pow(a * b)).toEqual(e2);
  });
}
