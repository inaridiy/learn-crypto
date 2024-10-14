import { BLS12_381_G1, BLS12_381_G2, CurvePoint, Polynomial } from "../paring";
import { FQ, curveOrder, fieldModulus } from "../paring/curve/bls12-381";

type Setup = { g1: CurvePoint; g2: CurvePoint }[];

const trustedSetup = (n: number, secret: bigint): Setup => {
  const setup: Setup = [];
  let s = 1n;
  for (let i = 0n; i < n; i++) {
    console.log(s);
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

const p2 = new Polynomial(
  [1n, 2n, 3n, 4n, 7n, 7n, 7n, 7n, 13n, 13n, 13n, 13n, 13n, 13n, 13n, 13n],
  curveOrder
);

{
  const input = 10n;
  const e = p2.eval(input);
  const setup = trustedSetup(p2.degree() + 1, input);
  const e1 = evalWithG1Setup(p2, setup);
  const e2 = BLS12_381_G1.mul(e.n);
  console.log(e1.x.eq(e2.x)); // ok true?
}

{
  const input = 1_927_409_816n;
  const e = p2.eval(input);
  const setup = trustedSetup(p2.degree() + 1, input);
  const e1 = evalWithG1Setup(p2, setup);
  const e2 = BLS12_381_G1.mul(e.n);
  console.log(e1.x.eq(e2.x)); // why false?
}
