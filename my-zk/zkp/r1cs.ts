import { FieldElement } from "../primitive/interface";

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
) => {};
