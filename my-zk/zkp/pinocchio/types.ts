import { ECPointCyclicGroupOnExtendedFF } from "../../primitive/elliptic-curve";
import { FiniteField, FiniteFieldElement } from "../../primitive/finite-field";
import { PolynomialOnFF } from "../../primitive/polynomial";
import { R1CSConstraints, StructuralWitness } from "../r1cs";

export type PinocchioPoint = {
  g1: ECPointCyclicGroupOnExtendedFF;
  g2: ECPointCyclicGroupOnExtendedFF;
};

export type PairingResult = {
  eq: (other: any) => boolean;
  mul: (other: any) => PairingResult;
};

export type PairingFunction = (
  g1: ECPointCyclicGroupOnExtendedFF,
  g2: ECPointCyclicGroupOnExtendedFF
) => PairingResult;

export type PinocchioGenerators = {
  g1: ECPointCyclicGroupOnExtendedFF;
  g2: ECPointCyclicGroupOnExtendedFF;
};

export type SampleFieldElement = (field: FiniteField) => FiniteFieldElement;

export type PinocchioConfig = {
  generators: PinocchioGenerators;
  pairing: PairingFunction;
  sampleFieldElement: SampleFieldElement;
};

export type PinocchioPointOps = {
  encode: (v: FiniteFieldElement) => PinocchioPoint;
  zero: () => PinocchioPoint;
  add: (a: PinocchioPoint, b: PinocchioPoint) => PinocchioPoint;
  scale: (point: PinocchioPoint, scalar: FiniteFieldElement) => PinocchioPoint;
};

export type PinocchioRuntime = {
  pairing: PairingFunction;
  sampleFieldElement: SampleFieldElement;
  pointOps: PinocchioPointOps;
};

/**
 * 線形結合チェック（α-Check）のためのCRS（共通参照文字列）部品。
 * q-PKE (q-Power Knowledge of Exponent) の仮定に基づき、
 * 証明者が多項式の線形結合を正しく計算したか検証するために用います。
 */
export type AlphaCrs = {
  gPolyAtS: PinocchioPoint[]; // { g^{u_i(s)} }
  gAlphaPolyAtS: PinocchioPoint[]; // { g^{\alpha \cdot u_i(s)} }
};

export type BetaBundles = PinocchioPoint[];

// QAP化された算術回路情報
export type PinocchioCircuit = {
  aps: PolynomialOnFF[]; // A行列由来の多項式群
  bps: PolynomialOnFF[]; // B行列由来の多項式群
  cps: PolynomialOnFF[]; // C行列由来の多項式群
  t: PolynomialOnFF; // ターゲット多項式
};

export type PinocchioKey = {
  gSs: PinocchioPoint[]; // sのべき乗のCRS {g, g^s, g^{s^2}, ...}
  crsApsMid: AlphaCrs; // 中間変数用 A多項式のαチェック用CRS
  crsBpsMid: AlphaCrs; // 中間変数用 B多項式のαチェック用CRS
  crsCpsMid: AlphaCrs; // 中間変数用 C多項式のαチェック用CRS
  betaBundles: BetaBundles; // βチェック用バンドルCRS
  gOne: PinocchioPoint; // g^1
  gAlphaA: PinocchioPoint; // g^{\alpha_A}
  gAlphaB: PinocchioPoint; // g^{\alpha_B}
  gAlphaC: PinocchioPoint; // g^{\alpha_C}
  gBeta: PinocchioPoint; // g^\beta
  gTPoly: PinocchioPoint; // g^{t(s)}
};

export type PinocchioMetadata = {
  field: FiniteField;
  inputCount: number;
  outputCount: number;
  privateInputCount: number;
  midStart: number; // witnessベクトルにおける秘密入力+中間変数の開始インデックス
  maxDegree: number;
  constraintCount: number;
};

// Trusted Setupの結果
export type PinocchioSetup = {
  circuit: PinocchioCircuit;
  key: PinocchioKey;
  metadata: PinocchioMetadata;
};

// Trusted Setupで生成される鍵（証明鍵・検証鍵の混合）
export type ZKPinocchioKey = {
  gSs: PinocchioPoint[]; // sのべき乗のCRS {g, g^s, g^{s^2}, ...}
  crsApsMid: AlphaCrs; // 中間変数用 A多項式のαチェック用CRS
  crsBpsMid: AlphaCrs; // 中間変数用 B多項式のαチェック用CRS
  crsCpsMid: AlphaCrs; // 中間変数用 C多項式のαチェック用CRS
  betaBundles: BetaBundles; // βチェック用バンドルCRS
  gOne: PinocchioPoint; // g^1
  gAlphaA: PinocchioPoint; // g^{\alpha_A}
  gAlphaB: PinocchioPoint; // g^{\alpha_B}
  gAlphaC: PinocchioPoint; // g^{\alpha_C}
  gAlphaAT: PinocchioPoint; // g^{\alpha_A t(s)}
  gAlphaBT: PinocchioPoint; // g^{\alpha_B t(s)}
  gAlphaCT: PinocchioPoint; // g^{\alpha_C t(s)}
  gBeta: PinocchioPoint; // g^\beta
  gBetaT: PinocchioPoint; // g^{\beta t(s)}
  gTPoly: PinocchioPoint; // g^{t(s)}
};

export type ZkPinocchioSetup = {
  circuit: PinocchioCircuit;
  key: ZKPinocchioKey;
  metadata: PinocchioMetadata;
};
/**
 * 証明オブジェクト π
 * 非常に簡潔（Succinct）で、定数個の群要素のみで構成されます。
 */
export type PinocchioProof = {
  gH: PinocchioPoint; // g^{h(s)}: 商多項式の評価値
  gAMid: PinocchioPoint; // g^{A_{mid}(s)}: 中間変数に関するAの線形結合
  gAlphaAMid: PinocchioPoint; // g^{\alpha_a A_{mid}(s)}: Aの線形結合チェック用
  gBMid: PinocchioPoint; // g^{B_{mid}(s)}
  gAlphaBMid: PinocchioPoint; // g^{\alpha_b B_{mid}(s)}
  gCMid: PinocchioPoint; // g^{C_{mid}(s)}
  gAlphaCMid: PinocchioPoint; // g^{\alpha_c C_{mid}(s)}
  gZ: PinocchioPoint; // g^Z: 係数一貫性（β）チェック用
};

export type ProveResult = {
  ioCoeffs: FiniteFieldElement[]; // 入出力（Public Input）
  proof: PinocchioProof; // 証明（Proof）
};

export type PinocchioInstance = {
  trustedSetup: (r1cs: R1CSConstraints<FiniteFieldElement>) => PinocchioSetup;
  prove: (setup: PinocchioSetup, witness: StructuralWitness<FiniteFieldElement>) => ProveResult;
  zkTrustedSetup: (r1cs: R1CSConstraints<FiniteFieldElement>) => ZkPinocchioSetup;
  zkProve: (setup: ZkPinocchioSetup, witness: StructuralWitness<FiniteFieldElement>) => ProveResult;
  verify: (setup: PinocchioSetup, proveResult: ProveResult) => boolean;
};
