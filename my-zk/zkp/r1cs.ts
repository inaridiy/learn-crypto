import { FieldElement } from "../primitive/interface";
import { FiniteField, FiniteFieldElement } from "../primitive/finite-field";

export type R1CSConstraints<F extends FieldElement<F, any>> = {
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
  const toFF = (value: number | bigint) =>
    field.from(typeof value === "bigint" ? value : BigInt(value));
  const row = (coeffs: Array<number | bigint>): FiniteFieldElement[] => coeffs.map(toFF);

  const r1cs: R1CSConstraints<FiniteFieldElement> = {
    a: [row([0, 1, 0, 0])],
    b: [row([0, 0, 1, 0])],
    c: [row([0, 0, 0, 1])],
  };

  const buildWitness = (
    x: number,
    y: number,
    z: number
  ): StructuralWitness<FiniteFieldElement> => ({
    one: field.one(),
    inputs: [toFF(x)],
    outputs: [toFF(y)],
    intermediates: [toFF(z)],
  });

  describe("isSatisfied", () => {
    it("returns true when an R1CS witness meets all constraints", () => {
      const witness = buildWitness(3, 4, 12);
      expect(isSatisfied(r1cs, witness)).toBe(true);
    });

    it("returns false when the witness violates a constraint", () => {
      const witness = buildWitness(3, 4, 10);
      expect(isSatisfied(r1cs, witness)).toBe(false);
    });

    it("throws when the witness does not start with the constant one", () => {
      const witness = buildWitness(2, 5, 10);
      witness.one = toFF(0);
      expect(() => isSatisfied(r1cs, witness)).toThrow("Witness must include the constant 1");
    });
  });
}
