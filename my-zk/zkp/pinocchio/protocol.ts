import { matrixToPolynomials } from "../../primitive/lagrange";
import { PolynomialFactory, PolynomialOnFF } from "../../primitive/polynomial";
import { FiniteFieldElement } from "../../primitive/finite-field";
import { ProveResult, PinocchioRuntime, PinocchioSetup } from "./types";
import { StructuralWitness, R1CSConstraints, toWitnessVector } from "../r1cs";
import { evalPolyWithCrs, linearCombinePolynomials, linearCombinationPoints } from "./operations";
import { buildAlphaCrs, buildBetaBundles, buildTargetPolynomial } from "./crs";

/*
 * Trusted Setup: トラステッドセットアップ
 * R1CSからQAPへの変換を行い、CRS（評価鍵・検証鍵）を生成します。
 * この関数は信頼できる第三者(TTP)によって一度だけ実行される想定です。
 */
export const createTrustedSetup =
  (runtime: PinocchioRuntime) =>
  (r1cs: R1CSConstraints<FiniteFieldElement>): PinocchioSetup => {
    const { pointOps, sampleFieldElement } = runtime;

    // Toxic Waste（有害な廃棄物）:
    // 以下のランダム値 s, α, β は、CRS生成後に完全に破棄されなければなりません。
    // これらが漏洩すると、任意の偽の証明を作成できてしまいます（健全性が崩壊する）。
    const [s, alphaA, alphaB, alphaC, beta] = Array.from({ length: 5 }, () =>
      sampleFieldElement(r1cs.structure)
    );

    const midStart = r1cs.index.output + 1; // 1, public inputs, public outputs の次から中間変数

    // FlatteningされたR1CSの行列A, B, Cを、ラグランジュ補間を用いて多項式の列に変換します（QAP化）。
    const [aps, bps, cps] = [r1cs.a, r1cs.b, r1cs.c].map((v) =>
      matrixToPolynomials(r1cs.structure, v)
    );

    // 中間変数(Mid)部分に対してのみ、健全性チェック用のCRSを生成します。
    // （入出力部分は検証者が自分で計算できるためCRSには含めません）
    const crsApsMid = buildAlphaCrs(pointOps, aps.slice(midStart), s, alphaA);
    const crsBpsMid = buildAlphaCrs(pointOps, bps.slice(midStart), s, alphaB);
    const crsCpsMid = buildAlphaCrs(pointOps, cps.slice(midStart), s, alphaC);

    // A, B, C の係数一貫性をチェックするためのCRSバンドルを生成します。
    const betaBundles = buildBetaBundles(
      pointOps,
      [aps.slice(midStart), bps.slice(midStart), cps.slice(midStart)],
      s,
      beta
    );

    // ターゲット多項式 t(x) を構築します。t(x) = Π(x - i)
    const pFactory = new PolynomialFactory(r1cs.structure);
    const t = buildTargetPolynomial(pFactory, r1cs.a.length);

    // 多項式の次数に合わせて s のべき乗 {s^1, s^2, ...} を用意します。
    const degreeCandidates = [...aps, ...bps, ...cps].map((p) => p.degree());
    const baseDegree = degreeCandidates.length > 0 ? Math.max(...degreeCandidates) : 0;
    const maxDegree = Math.max(baseDegree, t.degree());

    const ss = Array.from({ length: maxDegree }, (_, i) => s.pow(BigInt(i + 1)));
    const gSs = ss.map((value) => pointOps.encode(value));

    // 検証に必要な定数などをエンコード
    const gOne = pointOps.encode(r1cs.structure.one());
    const gAlphaA = pointOps.encode(alphaA);
    const gAlphaB = pointOps.encode(alphaB);
    const gAlphaC = pointOps.encode(alphaC);
    const gBeta = pointOps.encode(beta);
    const gTPoly = pointOps.encode(t.eval(s));

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

    const metadata = {
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

/**
 * 証明生成（Proving）
 * Witness（秘密の入力や中間変数の値）を用いて、計算が正しいことを証明する証拠を作成します。
 * ※この実装はVC(Verifiable Computation)としてのPinocchioであり、
 *  ランダムシフトによる完全なゼロ知識化（Zero-Knowledge）までは含まれていません。
 */
export const createProve =
  (runtime: PinocchioRuntime) =>
  (setup: PinocchioSetup, witness: StructuralWitness<FiniteFieldElement>): ProveResult => {
    const { circuit, key, metadata } = setup;
    const { pointOps } = runtime;
    const { aps, bps, cps, t } = circuit;
    const midStart = metadata.midStart;

    // Witnessをベクトル形式に変換 w = [1, inputs..., outputs..., intermediates...]
    const witnessVector = toWitnessVector(witness);
    const ioCoeffs = witnessVector.slice(0, midStart); // 公開部分（IO）
    const midCoeffs = witnessVector.slice(midStart); // 秘密部分（Mid）

    const pFactory = new PolynomialFactory(metadata.field);

    const combine = (polys: PolynomialOnFF[], coeffs: FiniteFieldElement[]) =>
      linearCombinePolynomials(pFactory, polys, coeffs);

    // A(x), B(x), C(x) をIO部分とMid部分に分けて計算します。
    const aIOPoly = combine(aps.slice(0, midStart), ioCoeffs);
    const bIOPoly = combine(bps.slice(0, midStart), ioCoeffs);
    const cIOPoly = combine(cps.slice(0, midStart), ioCoeffs);

    const aMidPoly = combine(aps.slice(midStart), midCoeffs);
    const bMidPoly = combine(bps.slice(midStart), midCoeffs);
    const cMidPoly = combine(cps.slice(midStart), midCoeffs);

    // IOとMidを結合して、完全な A(x), B(x), C(x) を得ます。
    const aPoly = aIOPoly.add(aMidPoly);
    const bPoly = bIOPoly.add(bMidPoly);
    const cPoly = cIOPoly.add(cMidPoly);

    // QAPの核心となる方程式: P(x) = A(x) * B(x) - C(x)
    const pPoly = aPoly.mul(bPoly).sub(cPoly);

    // 因数定理より、すべてのゲートが正しければ P(x) は t(x) で割り切れるはずです。
    // 商多項式 h(x) = P(x) / t(x) を求めます。
    const [hPoly, remainder] = pPoly.divmod(t);

    // 割り切れなければ、Witnessが間違っている（不正な計算をしている）ことになります。
    if (!remainder.isZero()) throw new Error("Constraint polynomial is not divisible by target");

    // (1) 商多項式 g^{h(s)} の計算
    // h(x)の係数とCRS {g^{s^i}} を使って計算します。
    const gH = evalPolyWithCrs(pointOps, hPoly, key.gOne, key.gSs);

    // (2) 中間変数部分の多項式 A_{mid}(s) とそのαチェック項の計算
    const gAMid = linearCombinationPoints(pointOps, key.crsApsMid.gPolyAtS, midCoeffs);
    const gAlphaAMid = linearCombinationPoints(pointOps, key.crsApsMid.gAlphaPolyAtS, midCoeffs);

    // (3) 中間変数部分の多項式 B_{mid}(s) とそのαチェック項の計算
    const gBMid = linearCombinationPoints(pointOps, key.crsBpsMid.gPolyAtS, midCoeffs);
    const gAlphaBMid = linearCombinationPoints(pointOps, key.crsBpsMid.gAlphaPolyAtS, midCoeffs);

    // (4) 中間変数部分の多項式 C_{mid}(s) とそのαチェック項の計算
    const gCMid = linearCombinationPoints(pointOps, key.crsCpsMid.gPolyAtS, midCoeffs);
    const gAlphaCMid = linearCombinationPoints(pointOps, key.crsCpsMid.gAlphaPolyAtS, midCoeffs);

    // (5) 係数一貫性チェック用 g^Z の計算 (β-Check)
    // 全ての中間変数に対して、A, B, Cで同じ係数 w_k が使われていることを証明するための項です。
    // Z = Σ_{k \in Mid} w_k * \beta(A_k(s) + B_k(s) + C_k(s))
    const gZ = linearCombinationPoints(pointOps, key.betaBundles, midCoeffs);

    // IO部分と証明（Proof）を返却
    return {
      ioCoeffs,
      proof: { gH, gAMid, gAlphaAMid, gBMid, gAlphaBMid, gCMid, gAlphaCMid, gZ },
    };
  };

/**
 * 検証（Verification）
 * 受け取った証明πと公開入出力(IO)を用いて、計算が正しかったかを検証します。
 * ペアリングという「魔法の写像」を使うことで、定数回の計算で検証可能です。
 */
export const createVerify =
  (runtime: PinocchioRuntime) =>
  (setup: PinocchioSetup, proveResult: ProveResult): boolean => {
    const { circuit, key, metadata } = setup;
    const { aps, bps, cps } = circuit;
    const { ioCoeffs, proof } = proveResult;
    const midStart = metadata.midStart;

    // 検証者はIO部分（公開情報）については自分で多項式を計算し、CRSで評価できます。
    // これにより、証明者が「どの入出力に対して計算したか」を固定します。
    const pFactory = new PolynomialFactory(metadata.field);
    const combine = (polys: PolynomialOnFF[], coeffs: FiniteFieldElement[]) =>
      linearCombinePolynomials(pFactory, polys, coeffs);

    const aIOPoly = combine(aps.slice(0, midStart), ioCoeffs);
    const bIOPoly = combine(bps.slice(0, midStart), ioCoeffs);
    const cIOPoly = combine(cps.slice(0, midStart), ioCoeffs);

    const gAIo = evalPolyWithCrs(runtime.pointOps, aIOPoly, key.gOne, key.gSs);
    const gBIo = evalPolyWithCrs(runtime.pointOps, bIOPoly, key.gOne, key.gSs);
    const gCIo = evalPolyWithCrs(runtime.pointOps, cIOPoly, key.gOne, key.gSs);

    // --------------------------------------------------------------------------
    // Check 1: ペアリングによるQAPチェック (Divisibility Check)
    // 証明された多項式が正しく t(s) で割り切れているか、つまり
    // A(s) * B(s) - C(s) = h(s) * t(s) が成り立つかを検証します。
    // 式: e(g^A, g^B) ?= e(g^t, g^h) * e(g^C, g^1)
    // ※ペアリングの双線形性 e(g^a, g^b) = e(g, g)^{ab} を利用
    // --------------------------------------------------------------------------

    // 全体の多項式評価値 = IO部分 + Mid部分(証明者から提示)
    const gA = runtime.pointOps.add(gAIo, proof.gAMid);
    const gB = runtime.pointOps.add(gBIo, proof.gBMid);
    const gC = runtime.pointOps.add(gCIo, proof.gCMid);

    const pair = runtime.pairing;
    const term1 = pair(gA.g1, gB.g2); // e(g^A, g^B)
    const term2 = pair(key.gTPoly.g1, proof.gH.g2); // e(g^t, g^h)
    const term3 = pair(gC.g1, key.gOne.g2); // e(g^C, g^1)

    // A*B = h*t + C  <=>  A*B - C = h*t
    const qapCheck = term1.eq(term2.mul(term3));

    // --------------------------------------------------------------------------
    // Check 2: 線形結合チェック (α-Check)
    // 証明者が勝手に g^{A_mid} を作ったのではなく、
    // CRSに含まれる多項式の線形結合として正しく構成したかを検証します。
    // 式: e(g^{\alpha A'}, g^1) ?= e(g^{A'}, g^\alpha)
    // --------------------------------------------------------------------------

    const alphaCheckA = pair(proof.gAlphaAMid.g1, key.gOne.g2).eq(
      pair(proof.gAMid.g1, key.gAlphaA.g2)
    );

    const alphaCheckB = pair(proof.gAlphaBMid.g1, key.gOne.g2).eq(
      pair(proof.gBMid.g1, key.gAlphaB.g2)
    );

    const alphaCheckC = pair(proof.gAlphaCMid.g1, key.gOne.g2).eq(
      pair(proof.gCMid.g1, key.gAlphaC.g2)
    );

    // --------------------------------------------------------------------------
    // Check 3: 係数一貫性チェック (β-Check)
    // A, B, C の線形結合において、同じ係数 w_i が使われているかを検証します。
    // 式: e(g^Z, g^1) ?= e(g^{A'} * g^{B'} * g^{C'}, g^\beta)
    // --------------------------------------------------------------------------

    const gSumMid = runtime.pointOps.add(
      runtime.pointOps.add(proof.gAMid, proof.gBMid),
      proof.gCMid
    );
    const betaCheck = pair(proof.gZ.g1, key.gOne.g2).eq(pair(gSumMid.g1, key.gBeta.g2));

    // 全てのチェックに通れば、検証成功です。
    return qapCheck && alphaCheckA && alphaCheckB && alphaCheckC && betaCheck;
  };
