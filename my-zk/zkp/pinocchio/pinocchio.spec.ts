import { describe, expect, it } from "vitest";
import { CURVE_ORDER } from "../../curves/bls12-381";
import { FiniteField, FiniteFieldElement } from "../../primitive/finite-field";
import { defaultPinocchioConfig } from "../pinocchio";
import { R1CSConstraints, StructuralWitness, toWitnessVector } from "../r1cs";
import { pinocchio } from ".";

const { trustedSetup, prove, verify } = pinocchio(defaultPinocchioConfig);

const field = new FiniteField(CURVE_ORDER);
const row = (coeffs: (bigint | number)[]): FiniteFieldElement[] =>
  coeffs.map((v) => field.from(BigInt(v)));
const matrix = (m: (bigint | number)[][]) => m.map((coeffs) => row(coeffs));

// テスト用のR1CSを定義する。
// 制約: (1*in_1) * (1*in_2) = (1*out_1)
const r1cs: R1CSConstraints<FiniteFieldElement> = {
  structure: field,
  index: { one: 0, input: 1, output: 2 },
  a: [row([0n, 1n, 0n, 0n])],
  b: [row([0n, 0n, 1n, 0n])],
  c: [row([0n, 0n, 0n, 1n])],
};

// 複雑な計算: out = x^3 + x + 5
const complexR1cs: R1CSConstraints<FiniteFieldElement> = {
  structure: field,
  index: { one: 0, input: 1, output: 2 },
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
  intermediates: bigint[] // intermediate
): StructuralWitness<FiniteFieldElement> => ({
  one: field.one(),
  inputs: [field.from(x)],
  outputs: [field.from(y)],
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
    const witness = buildWitness(3n, 4n, [12n]);
    const setup = trustedSetup(r1cs);
    const { ioCoeffs: _, proof } = prove(setup, witness);

    // 公開情報を改ざん (outputを4から6へ)
    // 正しい証明を持っていても、主張するOutputが異なれば検証に失敗するはず
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
    // 3^3 + 3 + 5 = 27 + 3 + 5 = 35
    const witness = buildWitness(3n, 35n, [9n, 27n, 30n]);
    const setup = trustedSetup(complexR1cs);
    const proof = prove(setup, witness);

    expect(verify(setup, proof)).toBe(true);
  });

  it("fails verification when public outputs differ with complex r1cs", () => {
    const witness = buildWitness(3n, 35n, [9n, 27n, 30n]);
    const setup = trustedSetup(complexR1cs);
    const { proof } = prove(setup, witness);

    const tamperedWitness = { ...witness, outputs: [field.from(6n)] }; // 不正な出力

    const tampered = {
      ioCoeffs: toWitnessVector(tamperedWitness).slice(0, setup.metadata.midStart),
      proof,
    };

    expect(verify(setup, tampered)).toBe(false);
  });
});
