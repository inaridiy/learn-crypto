import { BLS12_381_G1, BLS12_381_G2, CURVE_ORDER } from "./curves/bls12-381";
import { ECPointOnExtendedFF } from "./primitive/elliptic-curve";
import { ExtendedFiniteField } from "./primitive/extended-finite-field";
import { FiniteField, FiniteFieldElement } from "./primitive/finite-field";
import { pairing } from "./primitive/paring";
import { PolynomialFactory, PolynomialOnFF } from "./primitive/polynomial";

export const FF = new FiniteField(CURVE_ORDER);
export const POLY_FOR_KZG = new PolynomialFactory(FF);

export type TrustedSetup = {
  g1: ECPointOnExtendedFF;
  g2: ECPointOnExtendedFF;
}[];

export const makeTrustedSetup = (n: number, seed: bigint): TrustedSetup => {
  const setup: TrustedSetup = [];
  let s = 1n;
  for (let i = 0n; i < n; i++) {
    setup.push({
      g1: BLS12_381_G1.generator().scale(s % CURVE_ORDER),
      g2: BLS12_381_G2.generator().scale(s % CURVE_ORDER),
    });
    s = s * seed;
  }
  return setup;
};

const commitWithG1 = (p: PolynomialOnFF, setup: TrustedSetup) => {
  const { coeffs } = p;
  return coeffs
    .map((coeff, index) => setup[index].g1.scale(coeff.n))
    .reduce((acc, curr) => acc.add(curr));
};

const makeProof = (p: PolynomialOnFF, a: FiniteFieldElement, setup: TrustedSetup) => {
  const y = p.eval(a);
  const qx = p.sub(POLY_FOR_KZG.from([y])).quotient(POLY_FOR_KZG.from([-a.n, 1n]));
  const proof = commitWithG1(qx, setup);
  return { proof, y };
};

const verify = (
  proof: ECPointOnExtendedFF,
  a: FiniteFieldElement,
  y: FiniteFieldElement,
  commitment: ECPointOnExtendedFF,
  setup: TrustedSetup
) => {
  const left = pairing(proof, setup[1].g2.sub(BLS12_381_G2.generator().scale(a.n)));
  const right = pairing(
    commitment.sub(BLS12_381_G1.generator().scale(y.n)),
    BLS12_381_G2.generator()
  );
  return left.eq(right);
};

if (import.meta.vitest) {
  const { it, describe } = import.meta.vitest;
  describe("Kate Commitment", () => {
    it("should commit to a polynomial", () => {
      const p = POLY_FOR_KZG.from([1n, 2n, 3n, 4n]);

      const setup = makeTrustedSetup(p.degree() + 1, 114514n);
      const commitment = commitWithG1(p, setup);

      const a = FF.from(17n);
      const { proof, y } = makeProof(p, a, setup);
      const isValid = verify(proof, a, y, commitment, setup);
      expect(isValid).toBe(true);
    });
  });
}
