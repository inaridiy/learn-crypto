import { FiniteFieldElement } from "../../primitive/finite-field";
import { PolynomialFactory, PolynomialOnFF } from "../../primitive/polynomial";
import { AlphaCrs, BetaBundles, PinocchioPointOps } from "./types";

// 多項式の列に対し、線形結合チェック(α-Check)用のCRSペアを構築します。
export const buildAlphaCrs = (
  pointOps: PinocchioPointOps,
  polynomials: PolynomialOnFF[],
  s: FiniteFieldElement,
  alpha: FiniteFieldElement
): AlphaCrs => {
  // Schwartz-Zippelの補題に基づき、ランダムな点 s で評価します。
  const evals = polynomials.map((poly) => poly.eval(s)); //[u_1(s), u_2(s), ...]

  // 通常の項 g^{u(s)}
  const gPolyAtS = evals.map((val) => pointOps.encode(val));
  // α倍した項 g^{\alpha \cdot u(s)}
  const gAlphaPolyAtS = gPolyAtS.map((point) => pointOps.scale(point, alpha));

  return { gPolyAtS, gAlphaPolyAtS };
};

// 係数の一貫性チェック(β-Check)用のCRSバンドルを構築します。
// A, B, C の多項式で同じ係数が使われていることを強制するために、これらをβ倍して束ねます。
// K = g^{\beta(A_k(s) + B_k(s) + C_k(s))}
export const buildBetaBundles = (
  pointOps: PinocchioPointOps,
  [aps, bps, cps]: [PolynomialOnFF[], PolynomialOnFF[], PolynomialOnFF[]],
  s: FiniteFieldElement,
  beta: FiniteFieldElement
): BetaBundles => {
  const bundles = Array.from({ length: aps.length }, (_, i) => [aps[i], bps[i], cps[i]] as const)
    .map(([ap, bp, cp]) => ap.eval(s).add(bp.eval(s)).add(cp.eval(s)).mul(beta)) // \beta(a_i(s) + b_i(s) + c_i(s))
    .map((v) => pointOps.encode(v));

  return bundles;
};

/**
 * ターゲット多項式 t(x) = (x-r_1)(x-r_2)...(x-r_n) を構築する。
 * QAPにおいて、この多項式 t(x) で p(x) が割り切れることが、
 * 元のR1CSの全ゲート（制約）が満たされていることと同値になります。
 */
export const buildTargetPolynomial = (
  factory: PolynomialFactory<FiniteFieldElement, bigint | FiniteFieldElement>,
  constraintCount: number
) => {
  let target = factory.one();
  for (let i = 0; i < constraintCount; i++) {
    // 根(root)は制約のインデックスに対応させます（例: 1, 2, 3...）
    const root = factory.coeffField.from(BigInt(i + 1));
    const factor = factory.from([root.negate(), factory.coeffField.one()]); // (x - root)
    target = target.mul(factor);
  }
  return target;
};
