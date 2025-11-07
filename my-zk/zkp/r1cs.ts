import { FieldElement } from "../primitive/interface";
import { FiniteField, FiniteFieldElement } from "../primitive/finite-field";

export type R1CSConstraints<F extends FieldElement<F, any>> = {
  structure: F["structure"];
  index: { one: 0; input: number; output: number };
  a: F[][];
  b: F[][];
  c: F[][];
};

export type StructuralWitness<F extends FieldElement<F, any>> = {
  one: F;
  inputs: F[];
  outputs: F[];
  intermediates: F[];
};

export const toWitnessVector = <F extends FieldElement<F, any>>(witness: StructuralWitness<F>) => {
  return [witness.one, ...witness.inputs, ...witness.outputs, ...witness.intermediates];
};

export const isSatisfied = <F extends FieldElement<F, any>>(
  r1cs: R1CSConstraints<F>,
  witness: StructuralWitness<F>
) => {
  const { a: matrixA, b: matrixB, c: matrixC } = r1cs;
  const constraintCount = matrixA.length;

  if (constraintCount !== matrixB.length || constraintCount !== matrixC.length) {
    throw new Error("R1CS matrices must have the same number of constraints");
  }

  if (!witness.one.isOne()) {
    throw new Error("Witness must include the constant 1 as the first entry");
  }

  const witnessVector = toWitnessVector(witness);
  const width = witnessVector.length;
  const zero = witness.one.structure.zero();

  const dotProduct = (row: F[]) => {
    let acc = zero;
    for (let col = 0; col < width; col++) {
      const coeff = row[col];
      if (!coeff || coeff.isZero()) continue;
      const wire = witnessVector[col];
      if (!wire) {
        throw new Error(`Witness is missing value for wire ${col}`);
      }
      acc = acc.add(coeff.mul(wire));
    }
    return acc;
  };

  for (let i = 0; i < constraintCount; i++) {
    const left = dotProduct(matrixA[i]);
    const right = dotProduct(matrixB[i]);
    const out = dotProduct(matrixC[i]);
    if (!left.mul(right).eq(out)) return false;
  }

  return true;
};

if (import.meta.vitest) {
  const { describe, it, expect } = import.meta.vitest;
  const field = new FiniteField(17n);
  // const toFF = (value: number | bigint) =>
  //   field.from(typeof value === "bigint" ? value : BigInt(value));
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

  describe("isSatisfied", () => {
    it("returns true when an R1CS witness meets all constraints", () => {
      const witness = buildWitness(3n, 4n, 12n);
      expect(isSatisfied(r1cs, witness)).toBe(true);
    });

    it("returns false when the witness violates a constraint", () => {
      const witness = buildWitness(3n, 4n, 10n);
      expect(isSatisfied(r1cs, witness)).toBe(false);
    });

    it("throws when the witness does not start with the constant one", () => {
      const witness = buildWitness(2n, 5n, 10n);
      witness.one = field.from(0n);
      expect(() => isSatisfied(r1cs, witness)).toThrow("Witness must include the constant 1");
    });
  });
}
