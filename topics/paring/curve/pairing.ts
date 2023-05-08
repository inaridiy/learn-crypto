import { print } from "../../../utils/print";
import { CurvePoint, PointCoord } from "./EllipticCurve";
import { FQ, FQ12, BLS12_381_FQ12, fieldModulus, curveOrder } from "./bls12-381";

//魔法の数字だぜ！！何をしたいのか全く分からないぜ！！！
const ateLoopCount = 15132376222941642752n;
const logAteLoopCount = 62n;

export const embedFQ12 = (P: CurvePoint) => {
  return new CurvePoint(FQ12.from(P.x), FQ12.from(P.y), BLS12_381_FQ12);
};

export const twist = (P: CurvePoint) => {
  if (P.x.modPoly.degree() !== 2) throw new Error("なんかよくわからんけどだめらしいぜ");
  const { x, y } = P;
  const xc = [x.value.coefficients[0].sub(x.value.coefficients[1]), x.value.coefficients[1]];
  const yc = [y.value.coefficients[0].sub(y.value.coefficients[1]), y.value.coefficients[1]];

  const nx = FQ12.from([xc[0].n, 0n, 0n, 0n, 0n, 0n, xc[1].n, 0n, 0n, 0n, 0n, 0n, 0n]);
  const ny = FQ12.from([yc[0].n, 0n, 0n, 0n, 0n, 0n, yc[1].n, 0n, 0n, 0n, 0n, 0n, 0n]);

  const mx = nx.div([0n, 0n, 1n]);
  const my = ny.div([0n, 0n, 0n, 1n]);

  return new CurvePoint(mx, my, BLS12_381_FQ12);
};

export const lineFunction = (P1: CurvePoint, P2: CurvePoint, Q: PointCoord) => {
  if (P1.x.eq(0n) || P2.x.eq(0n) || Q.x.eq(0n)) return FQ12.one();

  const { x: x1, y: y1 } = P1;
  const { x: x2, y: y2 } = P2;
  const { x: xq, y: yq } = Q;

  if (!x1.eq(x2)) {
    const m = y2.sub(y1).div(x2.sub(x1));
    return xq.sub(x1).mul(m).add(y1).sub(yq);
  } else if (y1.eq(y2)) {
    const m = x1.pow(2n).mul(3n).div(y1.mul(2n));
    return xq.sub(x1).mul(m).add(y1).sub(yq);
  } else {
    return xq.sub(x1);
  }
};

export const millerLoop = (P: CurvePoint, Q: CurvePoint) => {
  if (P.x.eq(0n) || Q.x.eq(0n)) return FQ12.one();

  let R = Q;
  let f = FQ12.one();
  for (let i = logAteLoopCount; i >= 0n; i--) {
    // if (i === 62n) console.log(f.toString());

    f = f.pow(2n).mul(lineFunction(R, R, P));
    R = R.mul(2n);

    if (ateLoopCount & (1n << i)) {
      f = f.mul(lineFunction(R, Q, P));
      R = R.add(Q);
    }
  }

  return f.pow((fieldModulus ** 12n - 1n) / curveOrder);
};

export const pairing = (P: CurvePoint, Q: CurvePoint) => {
  return millerLoop(embedFQ12(P), twist(Q));
};
