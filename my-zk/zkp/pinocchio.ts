import { random } from "../../topics/random";
import { BLS12_381_G1, BLS12_381_G2, CURVE_ORDER, FQ12 } from "../curves/bls12-381";
import { FiniteField, FiniteFieldElement } from "../primitive/finite-field";
import { matrixToPolynomials } from "../primitive/lagrange";
import { PolynomialFactory, PolynomialOnFF } from "../primitive/polynomial";
import { pairing } from "../primitive/paring";
import { R1CSConstraints, StructuralWitness, isSatisfied, toWitnessVector } from "./r1cs";
import { ECPointCyclicGroupOnExtendedFF } from "../primitive/elliptic-curve";

/**
 * @fileoverview Pinocchioプロトコルの実装。
 * 参考記事: https://zenn.dev/qope/articles/f94b37ff2d9541
 * この実装はzk-SNARKの一種であるPinocchioプロトコルを提供する。
 * Pinocchioプロトコルは、ある計算が正しく行われたことを、計算の詳細（秘匿情報）を明かすことなく、
 * コンパクトな証明（proof）を用いて検証可能にする技術である。
 *
 * プロトコルは主に3つのステップから構成される：
 * 1. Trusted Setup（信頼できる第三者によるセットアップ）:
 *    計算回路に依存する公開パラメータ（評価鍵と検証鍵）を生成する。
 *    このステップは一度だけ行われ、生成されたパラメータは再利用可能。
 *    `trustedSetup`関数がこれに相当する。
 *
 * 2. Prove（証明の生成）:
 *    Prover（証明者）が、公開情報と秘匿情報（witness）を用いて、計算が正しく行われたことの証明（proof）を生成する。
 *    `prove`関数がこれに相当する。
 *
 * 3. Verify（証明の検証）:
 *    Verifier（検証者）が、公開情報、証明、検証鍵を用いて、証明が正しいかを検証する。
 *    `verify`関数がこれに相当する。
 */

// Pinocchioプロトコルでは、G1とG2という2つの異なる巡回群の点を用いる。
// これはペアリング演算 e(g1, g2) を行うためである。
export type PinocchioPoint = {
  g1: ECPointCyclicGroupOnExtendedFF;
  g2: ECPointCyclicGroupOnExtendedFF;
};

// G1, G2両方の群における単位元（ゼロ点）を返す。
const pointZero = (): PinocchioPoint => ({
  g1: BLS12_381_G1.zero(),
  g2: BLS12_381_G2.zero(),
});

// G1, G2両方の群で点の加算を行う。
const pointAdd = (a: PinocchioPoint, b: PinocchioPoint): PinocchioPoint => ({
  g1: a.g1.add(b.g1),
  g2: a.g2.add(b.g2),
});

// G1, G2両方の群で点のスカラー倍算を行う。
const pointScale = (point: PinocchioPoint, scalar: FiniteFieldElement): PinocchioPoint => ({
  g1: point.g1.scale(scalar.n),
  g2: point.g2.scale(scalar.n),
});

/**
 * 複数の点（PinocchioPoint）と係数（FiniteFieldElement）の線形結合を計算する。
 * acc = Σ (points[i] * coeffs[i])
 * @param points - 計算対象の点の配列。
 * @param coeffs - 各点に乗算する係数の配列。
 * @returns 線形結合の結果の点。
 */
const linearCombinePoints = (
  points: PinocchioPoint[],
  coeffs: FiniteFieldElement[]
): PinocchioPoint => {
  let acc = pointZero();
  for (let i = 0; i < points.length; i++) {
    const coeff = coeffs[i];
    if (!coeff || coeff.isZero()) continue;
    acc = pointAdd(acc, pointScale(points[i], coeff));
  }
  return acc;
};

/**
 * 複数の多項式と係数の線形結合を計算する。
 * acc = Σ (polys[i] * coeffs[i])
 * @param factory - 多項式を生成するためのファクトリ。
 * @param polys - 計算対象の多項式の配列。
 * @param coeffs - 各多項式に乗算する係数の配列。
 * @returns 線形結合の結果の多項式。
 */
const linearCombinePolynomials = (
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
 * ターゲット多項式 t(x) = (x-1)(x-2)...(x-n) を構築する。
 * この多項式は、R1CSの各制約が満たされるべき点（1, 2, ..., n）を根に持つ。
 * p(x) が t(x) で割り切れることは、全ての制約が満たされていることと同値である。
 * @param factory - 多項式を生成するためのファクトリ。
 * @param constraintCount - R1CSの制約の数。
 * @returns ターゲット多項式 t(x)。
 */
const buildTargetPolynomial = (
  factory: PolynomialFactory<FiniteFieldElement, bigint | FiniteFieldElement>,
  constraintCount: number
) => {
  let target = factory.one();
  for (let i = 0; i < constraintCount; i++) {
    const root = factory.coeffField.from(BigInt(i + 1));
    const factor = factory.from([root.negate(), factory.coeffField.one()]); // (x - root)
    target = target.mul(factor);
  }
  return target;
};

// Prover（証明者）が証明を生成するために必要な情報。Trusted Setupで生成される。
// E(s^i), E(α*s^i) など、秘密の点sにおける多項式の評価値（をエンコードしたもの）が含まれる。
export type PinocchioEvaluationKey = {
  ts: PinocchioPoint[]; // [E(s), E(s^2), ..., E(s^d)]
  tAlphaS: PinocchioPoint[]; // [E(αs), E(αs^2), ..., E(αs^d)]
  tasMid: PinocchioPoint[]; // [E(v_mid_k(s))], k > N
  tbsMid: PinocchioPoint[]; // [E(w_mid_k(s))], k > N
  tcsMid: PinocchioPoint[]; // [E(y_mid_k(s))], k > N
  tBetaAasMid: PinocchioPoint[]; // [E(βa*v_mid_k(s))]
  tBetaBbsMid: PinocchioPoint[]; // [E(βb*w_mid_k(s))]
  tBetaCcsMid: PinocchioPoint[]; // [E(βc*y_mid_k(s))]
  tBetaAtTarget: PinocchioPoint; // E(βa*t(s))
  tBetaBtTarget: PinocchioPoint; // E(βb*t(s))
  tBetaCtTarget: PinocchioPoint; // E(βc*t(s))
};

// Verifier（検証者）が証明を検証するために必要な情報。Trusted Setupで生成される。
// 公開情報（IO）に関する多項式の評価値や、検証に必要なその他のエンコードされた値が含まれる。
export type PinocchioVerificationKey = {
  tOne: PinocchioPoint; // E(1)
  tAlpha: PinocchioPoint; // E(α)
  tGamma: PinocchioPoint; // E(γ)
  tBetaAgamma: PinocchioPoint; // E(βa*γ)
  tBetaBgamma: PinocchioPoint; // E(βb*γ)
  tBetaCgamma: PinocchioPoint; // E(βc*γ)
  tTarget: PinocchioPoint; // E(t(s))
  tA0: PinocchioPoint; // E(v_0(s))
  tB0: PinocchioPoint; // E(w_0(s))
  tC0: PinocchioPoint; // E(y_0(s))
  tAsIO: PinocchioPoint[]; // [E(v_io_k(s))] for k=1..N
  tBsIO: PinocchioPoint[]; // [E(w_io_k(s))]
  tCsIO: PinocchioPoint[]; // [E(y_io_k(s))]
};

// プロトコル全体で共有されるメタデータ。
export type PinocchioMetadata = {
  field: FiniteField;
  inputCount: number;
  outputCount: number;
  midStart: number; // witnessベクトルにおける中間変数の開始インデックス
  maxDegree: number;
  constraintCount: number;
};

// Trusted Setupの結果をまとめた型。
export type PinocchioSetup = {
  evaluationKey: PinocchioEvaluationKey;
  verificationKey: PinocchioVerificationKey;
  metadata: PinocchioMetadata;
};

// Proverが生成する証明。Verifierに渡される。
// これらは全て楕円曲線上の点であり、witnessに関する情報を直接漏らすことはない。
export type PinocchioProof = {
  vMid: PinocchioPoint; // E(v_mid(s))
  wMid: PinocchioPoint; // E(w_mid(s))
  yMid: PinocchioPoint; // E(y_mid(s))
  h: PinocchioPoint; // E(h(s))
  alphaH: PinocchioPoint; // E(α*h(s))
  betaV: PinocchioPoint; // E(βa*v_mid(s))
  betaW: PinocchioPoint; // E(βb*w_mid(s))
  betaY: PinocchioPoint; // E(βc*y_mid(s))
};

/**
 * Trusted Setup: 評価鍵と検証鍵を生成する。
 * この関数は信頼できる第三者によって一度だけ実行される。
 * @param r1cs - 計算を表すR1CS (Rank-1 Constraint System)。
 * @returns 評価鍵、検証鍵、メタデータを含むセットアップオブジェクト。
 */
export const trustedSetup = (r1cs: R1CSConstraints<FiniteFieldElement>): PinocchioSetup => {
  // R1CSの行列A, B, Cを多項式v_k(x), w_k(x), y_k(x)に変換する。
  const [aps, bps, cps] = [r1cs.a, r1cs.b, r1cs.c].map((v) =>
    matrixToPolynomials(r1cs.structure, v)
  );

  const factory = new PolynomialFactory(r1cs.structure);
  // ターゲット多項式 t(x) を構築する。
  const targetPoly = buildTargetPolynomial(factory, r1cs.a.length);
  const degreeCandidates = [...aps, ...bps, ...cps].map((p) => p.degree());
  const baseDegree = degreeCandidates.length > 0 ? Math.max(...degreeCandidates) : 0;
  const maxDegree = Math.max(baseDegree, targetPoly.degree());

  // toxic waste: セットアップ後に破棄されなければならない秘密のパラメータ。
  // これらが漏洩すると、偽の証明を生成することが可能になる。
  const [s, alpha, betaA, betaB, betaC, gamma] = Array.from({ length: 6 }, () =>
    random(r1cs.structure.p)
  ).map((v) => r1cs.structure.from(v));

  // G1とG2の生成元を取得。
  const [g1, g2] = [BLS12_381_G1.generator(), BLS12_381_G2.generator()];
  // 有限体の要素 v を、G1とG2上の点 E(v) = (v*g1, v*g2) にエンコードする関数。
  const setupPoint = (v: FiniteFieldElement): PinocchioPoint => ({
    g1: g1.scale(v.n),
    g2: g2.scale(v.n),
  });

  // sのべき乗 s, s^2, ..., s^d を計算する。
  const ss = Array.from({ length: maxDegree }, (_, i) => s.pow(BigInt(i + 1)));

  // 検証鍵(VK)の一部を生成する。
  // これらは証明の検証に必要となる。
  const tOne = setupPoint(r1cs.structure.one()); // E(1)
  const tAlpha = setupPoint(alpha); // E(α)
  const tGamma = setupPoint(gamma); // E(γ)
  const tBetaAgamma = setupPoint(gamma.mul(betaA)); // E(γ*βa)
  const tBetaBgamma = setupPoint(gamma.mul(betaB)); // E(γ*βb)
  const tBetaCgamma = setupPoint(gamma.mul(betaC)); // E(γ*βc)

  // 評価鍵(EK)の一部を生成する。
  // これらはProverがh(s)を計算するために必要。
  const ts = ss.map((value) => setupPoint(value)); // [E(s), E(s^2), ...]
  const tAlphaS = ss.map((value) => setupPoint(alpha.mul(value))); // [E(αs), E(αs^2), ...]

  // 各多項式を秘密の点 s で評価し、その結果をエンコードする。
  const tas = aps.map((p) => setupPoint(p.eval(s))); // [E(v_0(s)), E(v_1(s)), ...]
  const tbs = bps.map((p) => setupPoint(p.eval(s))); // [E(w_0(s)), E(w_1(s)), ...]
  const tcs = cps.map((p) => setupPoint(p.eval(s))); // [E(y_0(s)), E(y_1(s)), ...]

  // Knowledge of Coefficient (KC) Assumptionのための値を生成する。
  // Proverが中間変数の多項式を正しく構築したことを検証するために使用する。
  const tBetaAas = aps.map((p) => setupPoint(p.eval(s).mul(betaA)));
  const tBetaBbs = bps.map((p) => setupPoint(p.eval(s).mul(betaB)));
  const tBetaCcs = cps.map((p) => setupPoint(p.eval(s).mul(betaC)));

  // ターゲット多項式を s で評価し、エンコードする。
  const targetEval = targetPoly.eval(s);
  const tTarget = setupPoint(targetEval); // E(t(s))
  const tBetaAtTarget = setupPoint(targetEval.mul(betaA));
  const tBetaBtTarget = setupPoint(targetEval.mul(betaB));
  const tBetaCtTarget = setupPoint(targetEval.mul(betaC));

  const midStart = r1cs.index.output + 1; // witnessにおける中間変数の開始位置

  // 評価鍵を構築する。
  // mid(中間変数)に関する値のみをProverに渡す。
  const ek: PinocchioEvaluationKey = {
    ts,
    tAlphaS,
    tasMid: tas.slice(midStart),
    tbsMid: tbs.slice(midStart),
    tcsMid: tcs.slice(midStart),
    tBetaAasMid: tBetaAas.slice(midStart),
    tBetaBbsMid: tBetaBbs.slice(midStart),
    tBetaCcsMid: tBetaCcs.slice(midStart),
    tBetaAtTarget,
    tBetaBtTarget,
    tBetaCtTarget,
  };

  // 検証鍵を構築する。
  // IO(公開入出力)に関する値と、検証に必要なその他の値をVerifierに渡す。
  const vk: PinocchioVerificationKey = {
    tOne,
    tAlpha,
    tGamma,
    tBetaAgamma,
    tBetaBgamma,
    tBetaCgamma,
    tTarget,
    tA0: tas[0], // E(v_0(s))
    tB0: tbs[0], // E(w_0(s))
    tC0: tcs[0], // E(y_0(s))
    tAsIO: tas.slice(1, midStart),
    tBsIO: tbs.slice(1, midStart),
    tCsIO: tcs.slice(1, midStart),
  };

  // メタデータを設定する。
  const metadata: PinocchioMetadata = {
    field: r1cs.structure,
    inputCount: r1cs.index.input - r1cs.index.one,
    outputCount: r1cs.index.output - r1cs.index.input,
    midStart,
    maxDegree,
    constraintCount: r1cs.a.length,
  };

  return { evaluationKey: ek, verificationKey: vk, metadata };
};

/**
 * 係数を使って、エンコードされたsのべき乗の級数評価を行う。
 * P(s) = c_0 + c_1*s + c_2*s^2 + ...
 * E(P(s)) = c_0*E(1) + c_1*E(s) + c_2*E(s^2) + ...
 *        = E(1)^c_0 * E(s)^c_1 * E(s^2)^c_2 * ...
 * @param coeffs - 多項式の係数 [c_0, c_1, ...]。
 * @param constantPoint - E(1) に相当する点。
 * @param series - [E(s), E(s^2), ...] の配列。
 * @returns E(P(s)) の計算結果。
 */
const buildSeriesEvaluation = (
  coeffs: FiniteFieldElement[],
  constantPoint: PinocchioPoint,
  series: PinocchioPoint[]
) => {
  let acc = pointZero();
  const constantCoeff = coeffs[0];
  if (constantCoeff && !constantCoeff.isZero()) {
    acc = pointAdd(acc, pointScale(constantPoint, constantCoeff));
  }
  for (let i = 1; i < coeffs.length; i++) {
    const coeff = coeffs[i];
    if (!coeff || coeff.isZero()) continue;
    const seriesIndex = i - 1;
    if (seriesIndex >= series.length) {
      throw new Error("Evaluation key does not cover polynomial degree");
    }
    acc = pointAdd(acc, pointScale(series[seriesIndex], coeff));
  }
  return acc;
};

// ゼロ点が含まれている場合にペアリング計算が失敗するのを防ぐラッパー関数。
// 点がゼロの場合、ペアリングの結果は乗法単位元(1)となる。
const safePairing = (p: ECPointCyclicGroupOnExtendedFF, q: ECPointCyclicGroupOnExtendedFF) => {
  if (p.isZero() || q.isZero()) return FQ12.one();
  return pairing(p, q);
};

/**
 * Prove: witnessを用いて証明を生成する。
 * @param r1cs - R1CS。
 * @param setup - Trusted Setupで生成されたパラメータ。
 * @param witness - 計算の入出力および中間値。
 * @returns Pinocchioの証明オブジェクト。
 */
export const prove = (
  r1cs: R1CSConstraints<FiniteFieldElement>,
  setup: PinocchioSetup,
  witness: StructuralWitness<FiniteFieldElement>
): PinocchioProof => {
  // まず、与えられたwitnessがR1CS制約を満たしているか確認する。
  if (!isSatisfied(r1cs, witness)) throw new Error("Witness does not satisfy R1CS constraints");

  const { evaluationKey: ek, verificationKey: vk, metadata } = setup;

  const factory = new PolynomialFactory(r1cs.structure);
  // R1CSを行列から多項式へ変換
  const [aps, bps, cps] = [r1cs.a, r1cs.b, r1cs.c].map((v) =>
    matrixToPolynomials(r1cs.structure, v)
  );
  const targetPoly = buildTargetPolynomial(factory, metadata.constraintCount);

  // witnessをベクトル形式に変換し、公開情報(IO)と中間変数(mid)に分割する。
  const witnessVector = toWitnessVector(witness);
  const midStart = metadata.midStart;
  const midCoeffs = witnessVector.slice(midStart);
  const ioCoeffs = witnessVector.slice(0, midStart);

  const combine = (polys: PolynomialOnFF[], coeffs: FiniteFieldElement[]) =>
    linearCombinePolynomials(factory, polys, coeffs);

  // v(x), w(x), y(x) をIO部分とmid部分に分けて計算する。
  const vIO = combine(aps.slice(0, midStart), ioCoeffs);
  const wIO = combine(bps.slice(0, midStart), ioCoeffs);
  const yIO = combine(cps.slice(0, midStart), ioCoeffs);

  const vMidPoly = combine(aps.slice(midStart), midCoeffs);
  const wMidPoly = combine(bps.slice(midStart), midCoeffs);
  const yMidPoly = combine(cps.slice(midStart), midCoeffs);

  // IOとmidを結合して、完全な v(x), w(x), y(x) を得る。
  const vPoly = vIO.add(vMidPoly);
  const wPoly = wIO.add(wMidPoly);
  const yPoly = yIO.add(yMidPoly);

  // p(x) = v(x) * w(x) - y(x) を計算する。
  const pPoly = vPoly.mul(wPoly).sub(yPoly);
  // p(x) を t(x) で割り、h(x) = p(x) / t(x) を求める。
  const [hPoly, remainder] = pPoly.divmod(targetPoly);
  // witnessが正しければ、余り(remainder)はゼロになるはず。
  if (!remainder.isZero()) throw new Error("Constraint polynomial is not divisible by target");

  const hCoeffs = hPoly.coeffs;

  // 評価鍵(EK)を用いて、証明の各要素を計算する。
  // これらはすべて秘密の点 s における評価値をエンコードしたものである。
  // Proverは s を知らないが、EKのおかげでこれらの値を計算できる。
  const vMidPoint = linearCombinePoints(ek.tasMid, midCoeffs); // E(v_mid(s))
  const wMidPoint = linearCombinePoints(ek.tbsMid, midCoeffs); // E(w_mid(s))
  const yMidPoint = linearCombinePoints(ek.tcsMid, midCoeffs); // E(y_mid(s))

  // KC検証用の値を計算
  const betaVPoint = linearCombinePoints(ek.tBetaAasMid, midCoeffs); // E(βa*v_mid(s))
  const betaWPoint = linearCombinePoints(ek.tBetaBbsMid, midCoeffs); // E(βb*w_mid(s))
  const betaYPoint = linearCombinePoints(ek.tBetaCcsMid, midCoeffs); // E(βc*y_mid(s))

  // h(s) と α*h(s) の評価値をエンコード
  const hPoint = buildSeriesEvaluation(hCoeffs, vk.tOne, ek.ts); // E(h(s))
  const alphaHPoint = buildSeriesEvaluation(hCoeffs, vk.tAlpha, ek.tAlphaS); // E(α*h(s))

  return {
    vMid: vMidPoint,
    wMid: wMidPoint,
    yMid: yMidPoint,
    h: hPoint,
    alphaH: alphaHPoint,
    betaV: betaVPoint,
    betaW: betaWPoint,
    betaY: betaYPoint,
  };
};

/**
 * Verify: 与えられた証明が正しいか検証する。
 * @param setup - Trusted Setupで生成されたパラメータ。
 * @param publicWitness - 公開情報（入力と出力）。
 * @param proof - Proverによって生成された証明。
 * @returns 検証が成功した場合はtrue、そうでなければfalse。
 */
export const verify = (
  setup: PinocchioSetup,
  publicWitness: StructuralWitness<FiniteFieldElement>,
  proof: PinocchioProof
) => {
  const { verificationKey: vk, metadata } = setup;

  // 公開witnessの形式が正しいか基本的なチェックを行う。
  if (!publicWitness.one.isOne()) return false;
  if (publicWitness.inputs.length !== metadata.inputCount) return false;
  if (publicWitness.outputs.length !== metadata.outputCount) return false;

  // 公開witnessをベクトル形式に変換する。
  const witnessIO = [publicWitness.one, ...publicWitness.inputs, ...publicWitness.outputs];
  const ioCoeffs = witnessIO.slice(1);

  // 検証鍵(VK)と公開情報(IO)を用いて、E(v_io(s)), E(w_io(s)), E(y_io(s))を計算する。
  const vIOPoint = linearCombinePoints(vk.tAsIO, ioCoeffs);
  const wIOPoint = linearCombinePoints(vk.tBsIO, ioCoeffs);
  const yIOPoint = linearCombinePoints(vk.tCsIO, ioCoeffs);

  // v, w, y の完全な評価値 E(v(s)), E(w(s)), E(y(s)) を復元する。
  // E(v(s)) = E(v_0(s)) + E(v_io(s)) + E(v_mid(s))
  const vAggregate = pointAdd(pointAdd(vk.tA0, vIOPoint), proof.vMid);
  const wAggregate = pointAdd(pointAdd(vk.tB0, wIOPoint), proof.wMid);
  const yAggregate = pointAdd(pointAdd(vk.tC0, yIOPoint), proof.yMid);

  // --- ここからペアリングを用いた一連の検証 ---

  // 1. h(x)が正しく構築されたかの検証 (Knowledge of Coefficient)
  // e(E(αh), E(1)) == e(E(h), E(α))
  const pairingAlphaLeft = safePairing(proof.alphaH.g1, vk.tOne.g2);
  const pairingAlphaRight = safePairing(proof.h.g1, vk.tAlpha.g2);
  if (!pairingAlphaLeft.eq(pairingAlphaRight)) return false;

  // 2. v_mid(x)が正しく構築されたかの検証 (Knowledge of Coefficient)
  // e(E(βa*v_mid), E(γ)) == e(E(v_mid), E(βa*γ))
  const pairingBetaVLeft = safePairing(proof.betaV.g1, vk.tGamma.g2);
  const pairingBetaVRight = safePairing(proof.vMid.g1, vk.tBetaAgamma.g2);
  if (!pairingBetaVLeft.eq(pairingBetaVRight)) return false;

  // 3. w_mid(x)が正しく構築されたかの検証 (Knowledge of Coefficient)
  // e(E(βb*w_mid), E(γ)) == e(E(w_mid), E(βb*γ))
  const pairingBetaWLeft = safePairing(proof.betaW.g1, vk.tGamma.g2);
  const pairingBetaWRight = safePairing(proof.wMid.g1, vk.tBetaBgamma.g2);
  if (!pairingBetaWLeft.eq(pairingBetaWRight)) return false;

  // 4. y_mid(x)が正しく構築されたかの検証 (Knowledge of Coefficient)
  // e(E(βc*y_mid), E(γ)) == e(E(y_mid), E(βc*γ))
  const pairingBetaYLeft = safePairing(proof.betaY.g1, vk.tGamma.g2);
  const pairingBetaYRight = safePairing(proof.yMid.g1, vk.tBetaCgamma.g2);
  if (!pairingBetaYLeft.eq(pairingBetaYRight)) return false;

  // 5. メインの制約 v(s)w(s) - y(s) = h(s)t(s) の検証
  // ペアリングの双線形性を用いて変形し、 e(v,w) / e(y,1) == e(h,t) をチェックする。
  const lhs = safePairing(vAggregate.g1, wAggregate.g2); // e(E(v(s)), E(w(s)))
  const rhs = safePairing(proof.h.g1, vk.tTarget.g2); // e(E(h(s)), E(t(s)))
  const denominator = safePairing(yAggregate.g1, vk.tOne.g2); // e(E(y(s)), E(1))
  if (!lhs.div(denominator).eq(rhs)) return false;

  // 全ての検証をパスした場合のみ、trueを返す。
  return true;
};

// --- 以下はテストコード ---
if (import.meta.vitest) {
  const { describe, it, expect } = import.meta.vitest;

  const field = new FiniteField(CURVE_ORDER);
  const row = (coeffs: bigint[]): FiniteFieldElement[] => coeffs.map((v) => field.from(v));

  // テスト用のR1CSを定義する。
  // 制約: (1*in_1) * (1*in_2) = (1*out_1)
  // 例: in_1=3, in_2=4 => out_1=12
  // witness vector: [one, in_1, in_2, out_1]
  // A = [0, 1, 0, 0]
  // B = [0, 0, 1, 0]
  // C = [0, 0, 0, 1]
  const r1cs: R1CSConstraints<FiniteFieldElement> = {
    structure: field,
    index: { one: 0, input: 1, output: 2 }, // このR1CSはinputが1つ、outputが1つ
    a: [row([0n, 1n, 0n, 0n])],
    b: [row([0n, 0n, 1n, 0n])], // コードの実装上、witness vectorは[one, input, output, intermediate]の順になる
    c: [row([0n, 0n, 0n, 1n])], // そのため、bは[0,0,1,0]ではなく、[0,0,0,1]に対応する中間変数を指すべきだが、このテストの構成上、x*y=zをx(input),y(output),z(intermediate)でマッピングしているため、このままで正しい。
  };

  const buildWitness = (
    x: bigint, // input
    y: bigint, // output: 本来は計算結果だが、このテストでは回路が単純なため、擬似的にinputとして扱う
    z: bigint // intermediate (x * y の結果)
  ): StructuralWitness<FiniteFieldElement> => ({
    one: field.one(),
    inputs: [field.from(x)],
    outputs: [field.from(y)], // 本来はoutputだが、テストの制約式c*w=zに対応させるため
    intermediates: [field.from(z)],
  });

  describe("pinocchio", () => {
    it("produces a proof that verifies", () => {
      // 3 * 4 = 12, 正しいwitness
      const witness = buildWitness(3n, 4n, 12n);
      const setup = trustedSetup(r1cs);
      const proof = prove(r1cs, setup, witness);

      expect(verify(setup, witness, proof)).toBe(true);
    });

    it("fails verification when public outputs differ", () => {
      // witnessは正しいが、検証時に異なる公開情報(outputs)を使う
      const witness = buildWitness(3n, 4n, 12n);
      const setup = trustedSetup(r1cs);
      const proof = prove(r1cs, setup, witness);

      // 公開情報を改ざん (outputを4から6へ)
      const tampered = {
        ...witness,
        outputs: [field.from(6n)],
      };

      expect(verify(setup, tampered, proof)).toBe(false);
    });

    it("rejects witnesses that do not satisfy the circuit", () => {
      const setup = trustedSetup(r1cs);
      // 2 * 5 != 9, 不正なwitness
      const badWitness = buildWitness(2n, 5n, 9n);
      // prove関数は内部でisSatisfiedを呼ぶため、ここでエラーがスローされるはず
      expect(() => prove(r1cs, setup, badWitness)).toThrow(
        "Witness does not satisfy R1CS constraints"
      );
    });
  });
}
