import { print } from "../../utils/print";
import { Polynomial, BLS12_381_G1, BLS12_381_G2, CurvePoint, pairing } from "../paring";
import { curveOrder } from "../paring/curve/bls12-381";
import { random } from "../random";

type Setup = { g1: CurvePoint; g2: CurvePoint }[];

const trustedSetup = (n: number, secret: bigint): Setup => {
  const setup: Setup = [];
  let s = 1n;
  for (let i = 0n; i < n; i++) {
    setup.push({
      g1: BLS12_381_G1.mul(s % curveOrder),
      g2: BLS12_381_G2.mul(s % curveOrder),
    });
    s = s * secret;
  }

  return setup;
};

const evalWithG1Setup = (p: Polynomial, setup: Setup) => {
  const { coefficients } = p;
  return coefficients
    .map((coeff, index) => setup[index].g1.mul(coeff.n))
    .reduce((acc, curr) => acc.add(curr));
};

const makeProof = (p: Polynomial, a: bigint, setup: Setup) => {
  const y = p.eval(a);
  const qx = p.sub(new Polynomial([y], y.p)).div(new Polynomial([-a, 1n], y.p));

  const proof = evalWithG1Setup(qx, setup);
  return { proof, y: y.n };
};

const verify = (proof: CurvePoint, a: bigint, y: bigint, commitment: CurvePoint, setup: Setup) => {
  const left = pairing(proof, setup[1].g2.sub(BLS12_381_G2.mul(a)));
  const right = pairing(commitment.sub(BLS12_381_G1.mul(y)), BLS12_381_G2);
  return left.eq(right);
};

if (import.meta.vitest) {
  it("KZG", () => {
    const p = new Polynomial([1n, 2n, 3n, 4n], curveOrder);

    const setup = trustedSetup(p.degree() + 1, 114514n);

    const commitment = evalWithG1Setup(p, setup);

    const a = 17n;
    const { proof, y } = makeProof(p, a, setup);

    const isValid = verify(proof, a, y, commitment, setup);

    expect(isValid).toBe(true);
  });

  it("KZG homomorphism", () => {
    const degree = 4;
    const p1 = new Polynomial([1n, 2n, 3n, 4n], curveOrder);
    const p2 = new Polynomial([5n, 6n, 7n, 8n], curveOrder);

    const setup = trustedSetup(degree + 1, random(curveOrder));

    const commitment1 = evalWithG1Setup(p1, setup);
    const commitment2 = evalWithG1Setup(p2, setup);

    const commitment1_2 = evalWithG1Setup(p1.add(p2), setup);

    const combinedCommitment = commitment1.add(commitment2);
    expect(commitment1_2.eq(combinedCommitment)).toBe(true);

    const a = 17n;
    const { proof: proof1, y: y1 } = makeProof(p1, a, setup);
    const { proof: proof2, y: y2 } = makeProof(p2, a, setup);
    const { proof: proof3, y: y3 } = makeProof(p1.add(p2), a, setup);
    const isValid1 = verify(proof1, a, y1, commitment1, setup);
    const isValid2 = verify(proof2, a, y2, commitment2, setup);
    const isValid3 = verify(proof3, a, y3, commitment1_2, setup);

    expect(isValid1).toBe(true);
    expect(isValid2).toBe(true);
    expect(isValid3).toBe(true);
    expect(y1 + y2).toBe(y3);
  });
}
