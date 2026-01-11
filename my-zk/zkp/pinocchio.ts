import { random } from "../../topics/random";
import { BLS12_381_G1, BLS12_381_G2, CURVE_ORDER, FQ12 } from "../curves/bls12-381";
import { FiniteField, FiniteFieldElement } from "../primitive/finite-field";
import { matrixToPolynomials } from "../primitive/lagrange";
import { PolynomialFactory, PolynomialOnFF } from "../primitive/polynomial";
import { pairing } from "../primitive/paring";
import { R1CSConstraints, StructuralWitness, isSatisfied, toWitnessVector } from "./r1cs";
import { ECPointCyclicGroupOnExtendedFF } from "../primitive/elliptic-curve";

const [g1, g2] = [BLS12_381_G1.generator(), BLS12_381_G2.generator()];

// Pinocchioプロトコルでは、G1とG2という2つの異なる巡回群の点を用いる。
// が、面倒くさいので両方ともセットで計算する
export type PinocchioPoint = {
  g1: ECPointCyclicGroupOnExtendedFF;
  g2: ECPointCyclicGroupOnExtendedFF;
};

const encode = (v: FiniteFieldElement): PinocchioPoint => ({
  g1: g1.scale(v.n),
  g2: g2.scale(v.n),
});

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

export type AlphaCrs = {
  gUss: PinocchioPoint[]; //
  gAlphaUss: PinocchioPoint[]; //
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
 * 点の配列と係数の配列から線形結合(Multi-Scalar Multiplication)を計算する
 * result = Σ (points[i] * coeffs[i])
 */
const linearCombinationPoints = (
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

//多項式の列に対し、AlphaCheck用のパラメータを構築する
const buildAlphaCrs = (
  polynomials: PolynomialOnFF[],
  s: FiniteFieldElement,
  alpha: FiniteFieldElement
): AlphaCrs => {
  const uss = polynomials.map((poly) => poly.eval(s)); //[u_1(s), u_2(s), ...]

  const gUss = uss.map((val) => encode(val)); // [g^{u_1(s)}, g^{u_2(s)}, ...]
  const gAlphaUss = gUss.map((point) => pointScale(point, alpha)); // [g^{\alpha \cdot u_1(s)}, g^{\alpha \cdot u_2(s)}, ...]

  return { gUss, gAlphaUss };
};

export type BetaBundles = PinocchioPoint[]; //

//多項式の列の組に対し、BetaCheck用のパラメータを構築する
const buildBetaBundles = (
  [aps, bps, cps]: [PolynomialOnFF[], PolynomialOnFF[], PolynomialOnFF[]],
  s: FiniteFieldElement,
  beta: FiniteFieldElement
): BetaBundles => {
  const bundles = Array.from({ length: aps.length }, (_, i) => [aps[i], bps[i], cps[i]] as const)
    .map(([ap, bp, cp]) => ap.eval(s).add(bp.eval(s)).add(cp.eval(s)).mul(beta)) // \beta(a_i(s) + b_i(s) + c_i(s))
    .map((v) => encode(v)); // g^{\beta(a_i(s) + b_i(s) + c_i(s))}

  return bundles;
};

/**
 * ターゲット多項式 t(x) = (x-1)(x-2)...(x-n) を構築する。
 * この多項式は、R1CSの各制約が満たされるべき点（1, 2, ..., n）を根に持つ。
 * p(x) が t(x) で割り切れることは、全ての制約が満たされていることと同値である。
 * @param factory - 多項式を生成するためのファクトリ。
 * @param constraintCount - R1CSの制約の数。
 * @returns ターゲット多項式 t(x)。
 */
const buildTargetPolynomial = (factory: PolynomialOnFF["structure"], constraintCount: number) => {
  let target = factory.one();
  for (let i = 0; i < constraintCount; i++) {
    const root = factory.coeffField.from(BigInt(i + 1));
    const factor = factory.from([root.negate(), factory.coeffField.one()]); // (x - root)
    target = target.mul(factor);
  }
  return target;
};

const evalPolyWithCrs = (poly: PolynomialOnFF, gOne: PinocchioPoint, gSs: PinocchioPoint[]) => {
  let acc = pointZero();
  const constantCoeff = poly.coeffs[0];
  if (constantCoeff && !constantCoeff.isZero()) {
    acc = pointAdd(acc, pointScale(gOne, constantCoeff));
  }
  for (let i = 1; i < poly.coeffs.length; i++) {
    const coeff = poly.coeffs[i];
    if (!coeff || coeff.isZero()) continue;
    const seriesIndex = i - 1;
    if (seriesIndex >= gSs.length) {
      throw new Error("Evaluation key does not cover polynomial degree");
    }
    acc = pointAdd(acc, pointScale(gSs[seriesIndex], coeff));
  }
  return acc;
};

export type PinocchioCircuit = {
  aps: PolynomialOnFF[];
  bps: PolynomialOnFF[];
  cps: PolynomialOnFF[];
  t: PolynomialOnFF;
};

export type PinocchioKey = {
  gSs: PinocchioPoint[];
  crsApsMid: AlphaCrs;
  crsBpsMid: AlphaCrs;
  crsCpsMid: AlphaCrs;
  betaBundles: BetaBundles;
  gOne: PinocchioPoint;
  gAlphaA: PinocchioPoint;
  gAlphaB: PinocchioPoint;
  gAlphaC: PinocchioPoint;
  gBeta: PinocchioPoint;
  gTPoly: PinocchioPoint;
};

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
  circuit: PinocchioCircuit;
  key: PinocchioKey;
  metadata: PinocchioMetadata;
};

/**
 * 証明オブジェクト π
 */
export type PinocchioProof = {
  gH: PinocchioPoint; // g^{h(s)}
  gAMid: PinocchioPoint; // g^{A_{mid}(s)}
  gAlphaAMid: PinocchioPoint; // g^{\alpha_a A_{mid}(s)}
  gBMid: PinocchioPoint; // g^{B_{mid}(s)}
  gAlphaBMid: PinocchioPoint; // g^{\alpha_b B_{mid}(s)}
  gCMid: PinocchioPoint; // g^{C_{mid}(s)}
  gAlphaCMid: PinocchioPoint; // g^{\alpha_c C_{mid}(s)}
  gZ: PinocchioPoint; // g^Z (係数一貫性チェック用)
};

export type ProveResult = {
  ioCoeffs: FiniteFieldElement[];
  proof: PinocchioProof;
};

/*
 * Trusted Setup: 評価鍵と検証鍵を生成する。
 * この関数は信頼できる第三者によって一度だけ実行される。
 * @param r1cs - 計算を表すR1CS (Rank-1 Constraint System)。
 * @returns 評価鍵、検証鍵、メタデータを含むセットアップオブジェクト。
 */
export const trustedSetup = (r1cs: R1CSConstraints<FiniteFieldElement>): PinocchioSetup => {
  // toxic waste: セットアップ後に破棄されなければならない秘密のパラメータ。
  // これらが漏洩すると、偽の証明を生成することが可能になる。
  const [s, alphaA, alphaB, alphaC, beta] = Array.from({ length: 5 }, () =>
    random(r1cs.structure.p)
  ).map((v) => r1cs.structure.from(v));

  const midStart = r1cs.index.output + 1;

  // R1CSの行列A, B, Cをラグランジュ補完で多項式の列APS, BPS, CPSに変換する
  const [aps, bps, cps] = [r1cs.a, r1cs.b, r1cs.c].map((v) =>
    matrixToPolynomials(r1cs.structure, v)
  );

  const crsApsMid = buildAlphaCrs(aps.slice(midStart), s, alphaA);
  const crsBpsMid = buildAlphaCrs(bps.slice(midStart), s, alphaB);
  const crsCpsMid = buildAlphaCrs(cps.slice(midStart), s, alphaC);

  const betaBundles = buildBetaBundles(
    [aps.slice(midStart), bps.slice(midStart), cps.slice(midStart)],
    s,
    beta
  );

  // ターゲット多項式 t(x) を構築する。
  const pFactory = new PolynomialFactory(r1cs.structure);
  const t = buildTargetPolynomial(pFactory, r1cs.a.length);
  const degreeCandidates = [...aps, ...bps, ...cps].map((p) => p.degree());
  const baseDegree = degreeCandidates.length > 0 ? Math.max(...degreeCandidates) : 0;
  const maxDegree = Math.max(baseDegree, t.degree());

  const ss = Array.from({ length: maxDegree }, (_, i) => s.pow(BigInt(i + 1)));
  const gSs = ss.map((value) => encode(value));

  const gOne = encode(r1cs.structure.one());
  const gAlphaA = encode(alphaA);
  const gAlphaB = encode(alphaB);
  const gAlphaC = encode(alphaC);
  const gBeta = encode(beta);
  const gTPoly = encode(t.eval(s));

  const pinocchioKey = {
    gSs,
    crsApsMid,
    crsBpsMid,
    crsCpsMid,
    betaBundles,
    gOne,
    gAlphaA,
    gAlphaB,
    gAlphaC,
    gBeta,
    gTPoly,
  };

  const metadata: PinocchioMetadata = {
    field: r1cs.structure,
    inputCount: r1cs.index.input - r1cs.index.one,
    outputCount: r1cs.index.output - r1cs.index.input,
    midStart,
    maxDegree,
    constraintCount: r1cs.a.length,
  };

  return {
    circuit: { aps, bps, cps, t },
    key: pinocchioKey,
    metadata,
  };
};

export const prove = (
  setup: PinocchioSetup,
  witness: StructuralWitness<FiniteFieldElement>
): ProveResult => {
  const { circuit, key, metadata } = setup;
  const { aps, bps, cps, t } = circuit;
  const midStart = metadata.midStart;

  const witnessVector = toWitnessVector(witness);
  const ioCoeffs = witnessVector.slice(0, midStart);
  const midCoeffs = witnessVector.slice(midStart);

  const pFactory = new PolynomialFactory(metadata.field);

  const combine = (polys: PolynomialOnFF[], coeffs: FiniteFieldElement[]) =>
    linearCombinePolynomials(pFactory, polys, coeffs);

  // A(x), B(x), C(x) をIO部分とmid部分に分けて計算する。
  const aIOPoly = combine(aps.slice(0, midStart), ioCoeffs);
  const bIOPoly = combine(bps.slice(0, midStart), ioCoeffs);
  const cIOPoly = combine(cps.slice(0, midStart), ioCoeffs);

  const aMidPoly = combine(aps.slice(midStart), midCoeffs);
  const bMidPoly = combine(bps.slice(midStart), midCoeffs);
  const cMidPoly = combine(cps.slice(midStart), midCoeffs);

  // IOとmidを結合して、完全な A(x), B(x), C(x) を得る。
  const aPoly = aIOPoly.add(aMidPoly);
  const bPoly = bIOPoly.add(bMidPoly);
  const cPoly = cIOPoly.add(cMidPoly);

  // p(x) = v(x) * w(x) - y(x) を計算する。
  const pPoly = aPoly.mul(bPoly).sub(cPoly);
  // p(x) を t(x) で割り、h(x) = p(x) / t(x) を求める。
  const [hPoly, remainder] = pPoly.divmod(t);
  // witnessが正しければ、余り(remainder)はゼロになるはず。
  if (!remainder.isZero()) throw new Error("Constraint polynomial is not divisible by target");

  // (1) g^{h(s)} の計算
  // h(x) = h_0 + h_1 x + ... なので、gOneとgSsを使って計算
  const gH = evalPolyWithCrs(hPoly, key.gOne, key.gSs);

  // (2) g^{A_{mid}(s)} と g^{\alpha_a A_{mid}(s)} の計算
  const gAMid = linearCombinationPoints(key.crsApsMid.gUss, midCoeffs);
  const gAlphaAMid = linearCombinationPoints(key.crsApsMid.gAlphaUss, midCoeffs);

  // (3) g^{B_{mid}(s)} と g^{\alpha_b B_{mid}(s)} の計算
  const gBMid = linearCombinationPoints(key.crsBpsMid.gUss, midCoeffs);
  const gAlphaBMid = linearCombinationPoints(key.crsBpsMid.gAlphaUss, midCoeffs);

  // (4) g^{C_{mid}(s)} と g^{\alpha_c C_{mid}(s)} の計算
  const gCMid = linearCombinationPoints(key.crsCpsMid.gUss, midCoeffs);
  const gAlphaCMid = linearCombinationPoints(key.crsCpsMid.gAlphaUss, midCoeffs);

  // (5) 係数一貫性チェック用 g^Z の計算
  // Z = Σ_{k \in Mid} c_k * \beta(A_k(s) + B_k(s) + C_k(s))
  const gZ = linearCombinationPoints(key.betaBundles, midCoeffs);

  return {
    ioCoeffs,
    proof: { gH, gAMid, gAlphaAMid, gBMid, gAlphaBMid, gCMid, gAlphaCMid, gZ },
  };
};

export const verify = (setup: PinocchioSetup, proveResult: ProveResult): boolean => {
  const { circuit, key, metadata } = setup;
  const { aps, bps, cps, t } = circuit;
  const { ioCoeffs, proof } = proveResult;
  const midStart = metadata.midStart;

  const pFactory = new PolynomialFactory(metadata.field);

  const combine = (polys: PolynomialOnFF[], coeffs: FiniteFieldElement[]) =>
    linearCombinePolynomials(pFactory, polys, coeffs);

  const aIOPoly = combine(aps.slice(0, midStart), ioCoeffs);
  const bIOPoly = combine(bps.slice(0, midStart), ioCoeffs);
  const cIOPoly = combine(cps.slice(0, midStart), ioCoeffs);

  const gAIo = evalPolyWithCrs(aIOPoly, key.gOne, key.gSs);
  const gBIo = evalPolyWithCrs(bIOPoly, key.gOne, key.gSs);
  const gCIo = evalPolyWithCrs(cIOPoly, key.gOne, key.gSs);

  // --------------------------------------------------------------------------
  // Check 1: ペアリングによるQAPチェック (Slide 100)
  // e(g^{A_io} * g^{A_mid}, g^{B_io} * g^{B_mid}) ?= e(g^{t(s)}, g^{h(s)}) * e(g^{C_io} * g^{C_mid}, g^1)
  // --------------------------------------------------------------------------

  const gA = pointAdd(gAIo, proof.gAMid);
  const gB = pointAdd(gBIo, proof.gBMid);
  const gC = pointAdd(gCIo, proof.gCMid);

  // QAPチェック
  const term1 = pairing(gA.g1, gB.g2);
  const term2 = pairing(key.gTPoly.g1, proof.gH.g2);
  const term3 = pairing(gC.g1, key.gOne.g2);

  const qapCheck = term1.eq(term2.mul(term3));

  // --------------------------------------------------------------------------
  // Check 2: 線形結合チェック (α-Check)
  // e(g^{\alpha A_{mid}}, g^1) ?= e(g^{A_{mid}}, g^\alpha)
  // --------------------------------------------------------------------------

  const alphaCheckA = pairing(proof.gAlphaAMid.g1, key.gOne.g2).eq(
    pairing(proof.gAMid.g1, key.gAlphaA.g2)
  );

  const alphaCheckB = pairing(proof.gAlphaBMid.g1, key.gOne.g2).eq(
    pairing(proof.gBMid.g1, key.gAlphaB.g2)
  );

  const alphaCheckC = pairing(proof.gAlphaCMid.g1, key.gOne.g2).eq(
    pairing(proof.gCMid.g1, key.gAlphaC.g2)
  );

  // --------------------------------------------------------------------------
  // Check 3: 係数一貫性チェック (β-Check)
  // e(g^Z, g^1) ?= e(g^{A_{mid}} * g^{B_{mid}} * g^{C_{mid}}, g^\beta)
  // --------------------------------------------------------------------------

  const gSumMid = pointAdd(pointAdd(proof.gAMid, proof.gBMid), proof.gCMid);
  const betaCheck = pairing(proof.gZ.g1, key.gOne.g2).eq(pairing(gSumMid.g1, key.gBeta.g2));

  return qapCheck && alphaCheckA && alphaCheckB && alphaCheckC && betaCheck;
};

// --- 以下はテストコード ---
if (import.meta.vitest) {
  const { describe, it, expect } = import.meta.vitest;

  const field = new FiniteField(CURVE_ORDER);
  const row = (coeffs: (bigint | number)[]): FiniteFieldElement[] =>
    coeffs.map((v) => field.from(BigInt(v)));
  const matrix = (m: (bigint | number)[][]) => m.map((coeffs) => row(coeffs));

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

  // out = x^3 + x + 5 の解を知っているか
  /**
   * x2 = x * x
   * x3 = x2 * x
   * x3_x = x3 + x
   * out = x3_x + 5
   */
  const complexR1cs: R1CSConstraints<FiniteFieldElement> = {
    structure: field,
    index: { one: 0, input: 1, output: 2 }, // このR1CSはinputが1つ、outputが1つ
    a: matrix([
      [0, 1, 0, 0, 0, 0],
      [0, 0, 0, 1, 0, 0],
      [0, 1, 0, 0, 1, 0],
      [5, 0, 0, 0, 0, 1],
    ]),
    b: matrix([
      [0, 1, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0],
      [1, 0, 0, 0, 0, 0],
      [1, 0, 0, 0, 0, 0],
    ]),
    c: matrix([
      [0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 1],
      [0, 0, 1, 0, 0, 0],
    ]),
  };

  const buildWitness = (
    x: bigint, // input
    y: bigint, // output
    intermediates: bigint[] // intermediate (x * y の結果)
  ): StructuralWitness<FiniteFieldElement> => ({
    one: field.one(),
    inputs: [field.from(x)],
    outputs: [field.from(y)], // 本来はoutputだが、テストの制約式c*w=zに対応させるため
    intermediates: row(intermediates),
  });

  describe("pinocchio", () => {
    it("produces a proof that verifies", () => {
      // 3 * 4 = 12, 正しいwitness
      const witness = buildWitness(3n, 4n, [12n]);
      const setup = trustedSetup(r1cs);
      const proveResult = prove(setup, witness);

      expect(verify(setup, proveResult)).toBe(true);
    });

    it("fails verification when public outputs differ", () => {
      // witnessは正しいが、検証時に異なる公開情報(outputs)を使う
      const witness = buildWitness(3n, 4n, [12n]);
      const setup = trustedSetup(r1cs);
      const { ioCoeffs: _, proof } = prove(setup, witness);

      // 公開情報を改ざん (outputを4から6へ)
      const tamperedWitness = {
        ...witness,
        outputs: [field.from(6n)],
      };

      const tampered = {
        ioCoeffs: toWitnessVector(tamperedWitness).slice(0, setup.metadata.midStart),
        proof,
      };

      expect(verify(setup, tampered)).toBe(false);
    });

    it("produces a proof that verifies with complex r1cs", () => {
      const witness = buildWitness(3n, 35n, [9n, 27n, 30n]);
      const setup = trustedSetup(complexR1cs);
      const proof = prove(setup, witness);

      expect(verify(setup, proof)).toBe(true);
    });

    it("fails verification when public outputs differ with complex r1cs", () => {
      const witness = buildWitness(3n, 35n, [9n, 27n, 30n]);

      const setup = trustedSetup(complexR1cs);
      const { proof } = prove(setup, witness);

      // 公開情報を改ざん (outputを4から6へ)
      const tamperedWitness = { ...witness, outputs: [field.from(6n)] };

      const tampered = {
        ioCoeffs: toWitnessVector(tamperedWitness).slice(0, setup.metadata.midStart),
        proof,
      };

      expect(verify(setup, tampered)).toBe(false);
    });
  });
}
