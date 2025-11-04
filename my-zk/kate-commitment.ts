import { BLS12_381_G1, BLS12_381_G2, CURVE_ORDER } from "./curves/bls12-381";
import { EllipticCurvePoint } from "./primitive/curve";
import { ExtendedFiniteField } from "./primitive/extended-finite-field";
import { FiniteField, FiniteFieldFactory } from "./primitive/finite-field";
import { pairing } from "./primitive/paring";
import { Polynomial, PolynomialFactory } from "./primitive/polynomial";

export const FF = new FiniteFieldFactory(CURVE_ORDER);
export const POLY_FOR_KZG = new PolynomialFactory(FF);

export type TrustedSetup = {
  g1: EllipticCurvePoint<ExtendedFiniteField>;
  g2: EllipticCurvePoint<ExtendedFiniteField>;
}[];

export const makeTrustedSetup = (n: number, seed: bigint): TrustedSetup => {
  const setup: TrustedSetup = [];
  let s = 1n;
  for (let i = 0n; i < n; i++) {
    setup.push({
      g1: BLS12_381_G1.mul(s % CURVE_ORDER),
      g2: BLS12_381_G2.mul(s % CURVE_ORDER),
    });
    s = s * seed;
  }
  return setup;
};

const commitWithG1 = (p: Polynomial, setup: TrustedSetup) => {
  const { coeffs } = p;
  return coeffs
    .map((coeff, index) => setup[index].g1.mul(coeff.n))
    .reduce((acc, curr) => acc.add(curr));
};

const makeProof = (p: Polynomial, a: FiniteField, setup: TrustedSetup) => {
  const y = p.eval(a);
  const qx = p.sub(POLY_FOR_KZG.from([y])).div(POLY_FOR_KZG.from([-a.n, 1n]));
  const proof = commitWithG1(qx, setup);
  return { proof, y };
};

const verify = (
  proof: EllipticCurvePoint<ExtendedFiniteField>,
  a: FiniteField,
  y: FiniteField,
  commitment: EllipticCurvePoint<ExtendedFiniteField>,
  setup: TrustedSetup
) => {
  const left = pairing(proof, setup[1].g2.sub(BLS12_381_G2.mul(a.n)));
  const right = pairing(commitment.sub(BLS12_381_G1.mul(y.n)), BLS12_381_G2);
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
