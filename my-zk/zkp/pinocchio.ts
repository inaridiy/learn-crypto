import { random } from "../../topics/random";
import { BLS12_381_G1, BLS12_381_G2, CURVE_ORDER, FQ12 } from "../curves/bls12-381";
import { FiniteField, FiniteFieldElement } from "../primitive/finite-field";
import { matrixToPolynomials } from "../primitive/lagrange";
import { PolynomialFactory, PolynomialOnFF } from "../primitive/polynomial";
import { pairing } from "../primitive/paring";
import { R1CSConstraints, StructuralWitness, isSatisfied, toWitnessVector } from "./r1cs";
import { ECPointCyclicGroupOnExtendedFF } from "../primitive/elliptic-curve";

export type PinocchioPoint = {
  g1: ECPointCyclicGroupOnExtendedFF;
  g2: ECPointCyclicGroupOnExtendedFF;
};

const pointZero = (): PinocchioPoint => ({
  g1: BLS12_381_G1.zero(),
  g2: BLS12_381_G2.zero(),
});

const pointAdd = (a: PinocchioPoint, b: PinocchioPoint): PinocchioPoint => ({
  g1: a.g1.add(b.g1),
  g2: a.g2.add(b.g2),
});

const pointScale = (point: PinocchioPoint, scalar: FiniteFieldElement): PinocchioPoint => ({
  g1: point.g1.scale(scalar.n),
  g2: point.g2.scale(scalar.n),
});

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

const buildTargetPolynomial = (
  factory: PolynomialFactory<FiniteFieldElement, bigint | FiniteFieldElement>,
  constraintCount: number
) => {
  let target = factory.one();
  for (let i = 0; i < constraintCount; i++) {
    const root = factory.coeffField.from(BigInt(i + 1));
    const factor = factory.from([root.negate(), factory.coeffField.one()]);
    target = target.mul(factor);
  }
  return target;
};

export type PinocchioEvaluationKey = {
  ts: PinocchioPoint[];
  tAlphaS: PinocchioPoint[];
  tasMid: PinocchioPoint[];
  tbsMid: PinocchioPoint[];
  tcsMid: PinocchioPoint[];
  tBetaAasMid: PinocchioPoint[];
  tBetaBbsMid: PinocchioPoint[];
  tBetaCcsMid: PinocchioPoint[];
  tBetaAtTarget: PinocchioPoint;
  tBetaBtTarget: PinocchioPoint;
  tBetaCtTarget: PinocchioPoint;
};

export type PinocchioVerificationKey = {
  tOne: PinocchioPoint;
  tAlpha: PinocchioPoint;
  tGamma: PinocchioPoint;
  tBetaAgamma: PinocchioPoint;
  tBetaBgamma: PinocchioPoint;
  tBetaCgamma: PinocchioPoint;
  tTarget: PinocchioPoint;
  tA0: PinocchioPoint;
  tB0: PinocchioPoint;
  tC0: PinocchioPoint;
  tAsIO: PinocchioPoint[];
  tBsIO: PinocchioPoint[];
  tCsIO: PinocchioPoint[];
};

export type PinocchioMetadata = {
  field: FiniteField;
  inputCount: number;
  outputCount: number;
  midStart: number;
  maxDegree: number;
  constraintCount: number;
};

export type PinocchioSetup = {
  evaluationKey: PinocchioEvaluationKey;
  verificationKey: PinocchioVerificationKey;
  metadata: PinocchioMetadata;
};

export type PinocchioProof = {
  vMid: PinocchioPoint;
  wMid: PinocchioPoint;
  yMid: PinocchioPoint;
  h: PinocchioPoint;
  alphaH: PinocchioPoint;
  betaV: PinocchioPoint;
  betaW: PinocchioPoint;
  betaY: PinocchioPoint;
};

export const trustedSetup = (r1cs: R1CSConstraints<FiniteFieldElement>): PinocchioSetup => {
  const [aps, bps, cps] = [r1cs.a, r1cs.b, r1cs.c].map((v) =>
    matrixToPolynomials(r1cs.structure, v)
  );

  const factory = new PolynomialFactory(r1cs.structure);
  const targetPoly = buildTargetPolynomial(factory, r1cs.a.length);
  const degreeCandidates = [...aps, ...bps, ...cps].map((p) => p.degree());
  const baseDegree = degreeCandidates.length > 0 ? Math.max(...degreeCandidates) : 0;
  const maxDegree = Math.max(baseDegree, targetPoly.degree());

  const [s, alpha, betaA, betaB, betaC, gamma] = Array.from({ length: 6 }, () =>
    random(r1cs.structure.p)
  ).map((v) => r1cs.structure.from(v));

  const [g1, g2] = [BLS12_381_G1.generator(), BLS12_381_G2.generator()];
  const setupPoint = (v: FiniteFieldElement): PinocchioPoint => ({
    g1: g1.scale(v.n),
    g2: g2.scale(v.n),
  });

  const ss = Array.from({ length: maxDegree }, (_, i) => s.pow(BigInt(i + 1)));

  const tOne = setupPoint(r1cs.structure.one());
  const tAlpha = setupPoint(alpha);
  const tGamma = setupPoint(gamma);
  const tBetaAgamma = setupPoint(gamma.mul(betaA));
  const tBetaBgamma = setupPoint(gamma.mul(betaB));
  const tBetaCgamma = setupPoint(gamma.mul(betaC));

  const ts = ss.map((value) => setupPoint(value));
  const tAlphaS = ss.map((value) => setupPoint(alpha.mul(value)));

  const tas = aps.map((p) => setupPoint(p.eval(s)));
  const tbs = bps.map((p) => setupPoint(p.eval(s)));
  const tcs = cps.map((p) => setupPoint(p.eval(s)));

  const tBetaAas = aps.map((p) => setupPoint(p.eval(s).mul(betaA)));
  const tBetaBbs = bps.map((p) => setupPoint(p.eval(s).mul(betaB)));
  const tBetaCcs = cps.map((p) => setupPoint(p.eval(s).mul(betaC)));

  const targetEval = targetPoly.eval(s);
  const tTarget = setupPoint(targetEval);
  const tBetaAtTarget = setupPoint(targetEval.mul(betaA));
  const tBetaBtTarget = setupPoint(targetEval.mul(betaB));
  const tBetaCtTarget = setupPoint(targetEval.mul(betaC));

  const midStart = r1cs.index.output + 1;
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

  const vk: PinocchioVerificationKey = {
    tOne,
    tAlpha,
    tGamma,
    tBetaAgamma,
    tBetaBgamma,
    tBetaCgamma,
    tTarget,
    tA0: tas[0],
    tB0: tbs[0],
    tC0: tcs[0],
    tAsIO: tas.slice(1, midStart),
    tBsIO: tbs.slice(1, midStart),
    tCsIO: tcs.slice(1, midStart),
  };

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

const safePairing = (p: ECPointCyclicGroupOnExtendedFF, q: ECPointCyclicGroupOnExtendedFF) => {
  if (p.isZero() || q.isZero()) return FQ12.one();
  return pairing(p, q);
};

export const prove = (
  r1cs: R1CSConstraints<FiniteFieldElement>,
  setup: PinocchioSetup,
  witness: StructuralWitness<FiniteFieldElement>
): PinocchioProof => {
  if (!isSatisfied(r1cs, witness)) throw new Error("Witness does not satisfy R1CS constraints");

  const { evaluationKey: ek, verificationKey: vk, metadata } = setup;

  const factory = new PolynomialFactory(r1cs.structure);
  const [aps, bps, cps] = [r1cs.a, r1cs.b, r1cs.c].map((v) =>
    matrixToPolynomials(r1cs.structure, v)
  );
  const targetPoly = buildTargetPolynomial(factory, metadata.constraintCount);

  const witnessVector = toWitnessVector(witness);
  const midStart = metadata.midStart;
  const midCoeffs = witnessVector.slice(midStart);
  const ioCoeffs = witnessVector.slice(0, midStart);

  const combine = (polys: PolynomialOnFF[], coeffs: FiniteFieldElement[]) =>
    linearCombinePolynomials(factory, polys, coeffs);

  const vIO = combine(aps.slice(0, midStart), ioCoeffs);
  const wIO = combine(bps.slice(0, midStart), ioCoeffs);
  const yIO = combine(cps.slice(0, midStart), ioCoeffs);

  const vMidPoly = combine(aps.slice(midStart), midCoeffs);
  const wMidPoly = combine(bps.slice(midStart), midCoeffs);
  const yMidPoly = combine(cps.slice(midStart), midCoeffs);

  const vPoly = vIO.add(vMidPoly);
  const wPoly = wIO.add(wMidPoly);
  const yPoly = yIO.add(yMidPoly);

  const pPoly = vPoly.mul(wPoly).sub(yPoly);
  const [hPoly, remainder] = pPoly.divmod(targetPoly);
  if (!remainder.isZero()) throw new Error("Constraint polynomial is not divisible by target");

  const hCoeffs = hPoly.coeffs;

  const vMidPoint = linearCombinePoints(ek.tasMid, midCoeffs);
  const wMidPoint = linearCombinePoints(ek.tbsMid, midCoeffs);
  const yMidPoint = linearCombinePoints(ek.tcsMid, midCoeffs);

  const betaVPoint = linearCombinePoints(ek.tBetaAasMid, midCoeffs);
  const betaWPoint = linearCombinePoints(ek.tBetaBbsMid, midCoeffs);
  const betaYPoint = linearCombinePoints(ek.tBetaCcsMid, midCoeffs);

  const hPoint = buildSeriesEvaluation(hCoeffs, vk.tOne, ek.ts);
  const alphaHPoint = buildSeriesEvaluation(hCoeffs, vk.tAlpha, ek.tAlphaS);

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

export const verify = (
  setup: PinocchioSetup,
  publicWitness: StructuralWitness<FiniteFieldElement>,
  proof: PinocchioProof
) => {
  const { verificationKey: vk, metadata } = setup;

  if (!publicWitness.one.isOne()) return false;
  if (publicWitness.inputs.length !== metadata.inputCount) return false;
  if (publicWitness.outputs.length !== metadata.outputCount) return false;

  const witnessIO = [publicWitness.one, ...publicWitness.inputs, ...publicWitness.outputs];
  const ioCoeffs = witnessIO.slice(1);

  const vIOPoint = linearCombinePoints(vk.tAsIO, ioCoeffs);
  const wIOPoint = linearCombinePoints(vk.tBsIO, ioCoeffs);
  const yIOPoint = linearCombinePoints(vk.tCsIO, ioCoeffs);

  const vAggregate = pointAdd(pointAdd(vk.tA0, vIOPoint), proof.vMid);
  const wAggregate = pointAdd(pointAdd(vk.tB0, wIOPoint), proof.wMid);
  const yAggregate = pointAdd(pointAdd(vk.tC0, yIOPoint), proof.yMid);

  const pairingAlphaLeft = safePairing(proof.alphaH.g1, vk.tOne.g2);
  const pairingAlphaRight = safePairing(proof.h.g1, vk.tAlpha.g2);
  if (!pairingAlphaLeft.eq(pairingAlphaRight)) return false;

  const pairingBetaVLeft = safePairing(proof.betaV.g1, vk.tGamma.g2);
  const pairingBetaVRight = safePairing(proof.vMid.g1, vk.tBetaAgamma.g2);
  if (!pairingBetaVLeft.eq(pairingBetaVRight)) return false;

  const pairingBetaWLeft = safePairing(proof.betaW.g1, vk.tGamma.g2);
  const pairingBetaWRight = safePairing(proof.wMid.g1, vk.tBetaBgamma.g2);
  if (!pairingBetaWLeft.eq(pairingBetaWRight)) return false;

  const pairingBetaYLeft = safePairing(proof.betaY.g1, vk.tGamma.g2);
  const pairingBetaYRight = safePairing(proof.yMid.g1, vk.tBetaCgamma.g2);
  if (!pairingBetaYLeft.eq(pairingBetaYRight)) return false;

  const lhs = safePairing(vAggregate.g1, wAggregate.g2);
  const rhs = safePairing(proof.h.g1, vk.tTarget.g2);
  const denominator = safePairing(yAggregate.g1, vk.tOne.g2);
  if (!lhs.div(denominator).eq(rhs)) return false;

  return true;
};

if (import.meta.vitest) {
  const { describe, it, expect } = import.meta.vitest;

  const field = new FiniteField(CURVE_ORDER);
  const row = (coeffs: bigint[]): FiniteFieldElement[] => coeffs.map((v) => field.from(v));

  const r1cs: R1CSConstraints<FiniteFieldElement> = {
    structure: field,
    index: { one: 0, input: 1, output: 2 },
    a: [row([0n, 1n, 0n, 0n])],
    b: [row([0n, 0n, 1n, 0n])],
    c: [row([0n, 0n, 0n, 1n])],
  };

  const buildWitness = (
    x: bigint,
    y: bigint,
    z: bigint
  ): StructuralWitness<FiniteFieldElement> => ({
    one: field.one(),
    inputs: [field.from(x)],
    outputs: [field.from(y)],
    intermediates: [field.from(z)],
  });

  describe("pinocchio", () => {
    it("produces a proof that verifies", () => {
      const witness = buildWitness(3n, 4n, 12n);
      const setup = trustedSetup(r1cs);
      const proof = prove(r1cs, setup, witness);

      expect(verify(setup, witness, proof)).toBe(true);
    });

    it("fails verification when public outputs differ", () => {
      const witness = buildWitness(3n, 4n, 12n);
      const setup = trustedSetup(r1cs);
      const proof = prove(r1cs, setup, witness);

      const tampered = {
        ...witness,
        outputs: [field.from(6n)],
      };

      expect(verify(setup, tampered, proof)).toBe(false);
    });

    it("rejects witnesses that do not satisfy the circuit", () => {
      const setup = trustedSetup(r1cs);
      const badWitness = buildWitness(2n, 5n, 9n);
      expect(() => prove(r1cs, setup, badWitness)).toThrow(
        "Witness does not satisfy R1CS constraints"
      );
    });
  });
}
