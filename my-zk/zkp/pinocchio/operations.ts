import { FiniteFieldElement } from "../../primitive/finite-field";
import { PolynomialFactory, PolynomialOnFF } from "../../primitive/polynomial";
import { PinocchioGenerators, PinocchioPoint, PinocchioPointOps } from "./types";

export const createPointOps = (generators: PinocchioGenerators): PinocchioPointOps => {
  const encode = (v: FiniteFieldElement) => ({
    g1: generators.g1.scale(v.n),
    g2: generators.g2.scale(v.n),
  });

  const zero = () => ({
    g1: generators.g1.structure.zero(),
    g2: generators.g2.structure.zero(),
  });

  const add = (a: PinocchioPoint, b: PinocchioPoint) => ({
    g1: a.g1.add(b.g1),
    g2: a.g2.add(b.g2),
  });

  const scale = (point: PinocchioPoint, scalar: FiniteFieldElement) => ({
    g1: point.g1.scale(scalar.n),
    g2: point.g2.scale(scalar.n),
  });

  return { encode, zero, add, scale };
};

/**
 * 複数の多項式と係数の線形結合を計算する（多項式同士の足し合わせ）。
 * acc = Σ (coeffs[i] * polys[i])
 */
export const linearCombinePolynomials = (
  factory: PolynomialFactory<FiniteFieldElement, bigint | FiniteFieldElement>,
  polys: PolynomialOnFF[],
  coeffs: FiniteFieldElement[]
) => {
  let acc = factory.zero();
  for (let i = 0; i < polys.length; i++) {
    const coeff = coeffs[i];
    if (!coeff || coeff.isZero()) continue;
    acc = acc.add(polys[i].scale(coeff.n));
  }
  return acc;
};

/**
 * 群上の点と係数の線形結合(Multi-Scalar Multiplication)を計算する。
 * 証明生成時に、CRSの点に対して係数を掛けて足し合わせる処理（準同型性を利用）です。
 * result = Σ (coeffs[i] * points[i])
 */
export const linearCombinationPoints = (
  pointOps: PinocchioPointOps,
  points: PinocchioPoint[],
  coeffs: FiniteFieldElement[]
): PinocchioPoint => {
  let acc = pointOps.zero();
  for (let i = 0; i < points.length; i++) {
    const coeff = coeffs[i];
    if (!coeff || coeff.isZero()) continue;
    acc = pointOps.add(acc, pointOps.scale(points[i], coeff));
  }
  return acc;
};

/**
 * CRSに含まれる {g^{s^0}, g^{s^1}, ...} を用いて、
 * 多項式 poly(x) を点 s で評価した値を指数部に持つ点 g^{poly(s)} を計算します。
 * 証明者は s を知らずに計算できます。
 */
export const evalPolyWithCrs = (
  pointOps: PinocchioPointOps,
  poly: PolynomialOnFF,
  gOne: PinocchioPoint,
  gSs: PinocchioPoint[]
) => {
  let acc = pointOps.zero();
  const constantCoeff = poly.coeffs[0];
  if (constantCoeff && !constantCoeff.isZero()) {
    acc = pointOps.add(acc, pointOps.scale(gOne, constantCoeff));
  }
  for (let i = 1; i < poly.coeffs.length; i++) {
    const coeff = poly.coeffs[i];
    if (!coeff || coeff.isZero()) continue;
    const seriesIndex = i - 1;
    if (seriesIndex >= gSs.length) {
      throw new Error("Evaluation key does not cover polynomial degree");
    }
    acc = pointOps.add(acc, pointOps.scale(gSs[seriesIndex], coeff));
  }
  return acc;
};
