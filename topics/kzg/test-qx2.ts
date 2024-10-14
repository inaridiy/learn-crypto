import { ProjPointType } from "@noble/curves/abstract/weierstrass";
import { bls12_381 } from "@noble/curves/bls12-381";
import { Polynomial } from "../paring";

bls12_381.G1.ProjectivePoint.BASE.multiply;

type Setup = { g1: ProjPointType<bigint> }[];

const trustedSetup = (n: number, secret: bigint): Setup => {
  const setup: Setup = [];
  let s = 1n;
  for (let i = 0n; i < n; i++) {
    setup.push({
      g1: bls12_381.G1.ProjectivePoint.BASE.multiply(s % bls12_381.fields.Fr.ORDER),
    });
    s = s * secret;
  }

  return setup;
};

const evalWithG1Setup = (p: Polynomial, setup: Setup) => {
  const { coefficients } = p;
  return coefficients
    .map((coeff, index) => setup[index].g1.multiply(coeff.n))
    .reduce((acc, curr) => acc.add(curr));
};

const p2 = new Polynomial(
  [1n, 2n, 3n, 4n, 7n, 7n, 7n, 7n, 13n, 13n, 13n, 13n, 13n, 13n, 13n, 13n],
  bls12_381.fields.Fr.ORDER
);

{
  const input = 10n;
  const e = p2.eval(input);
  const setup = trustedSetup(p2.degree() + 1, input);
  const e1 = evalWithG1Setup(p2, setup);
  const e2 = bls12_381.G1.ProjectivePoint.BASE.multiply(e.n);
  console.log(e1.equals(e2));
}

{
  const input = 1_927_409_816n;
  const e = p2.eval(input);
  const setup = trustedSetup(p2.degree() + 1, input);
  const e1 = evalWithG1Setup(p2, setup);
  const e2 = bls12_381.G1.ProjectivePoint.BASE.multiply(e.n);
  console.log(e1.equals(e2));
}
