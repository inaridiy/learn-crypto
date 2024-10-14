import { ProjPointType } from "@noble/curves/abstract/weierstrass";
import { bls12_381 } from "@noble/curves/bls12-381";
import { BLS12_381_G1, fieldModulus } from "../paring/curve/bls12-381";
import { CurvePoint, Polynomial } from "../paring";

bls12_381.G1.ProjectivePoint.BASE.multiply;

type Setup = { s: bigint; g1: ProjPointType<bigint> }[];

const trustedSetup = (n: number, secret: bigint): Setup => {
  const setup: Setup = [];
  let s = 1n;
  for (let i = 0n; i < n; i++) {
    setup.push({
      s: s % bls12_381.fields.Fr.ORDER,
      g1: bls12_381.G1.ProjectivePoint.BASE.multiply(s % bls12_381.fields.Fr.ORDER),
    });
    s = s * secret;
  }

  return setup;
};

type Setup2 = { s: bigint; g1: CurvePoint }[];

const mytrustedSetup = (n: number, secret: bigint): Setup2 => {
  const setup: Setup2 = [];
  let s = 1n;
  for (let i = 0n; i < n; i++) {
    setup.push({
      s: s % bls12_381.fields.Fr.ORDER,
      g1: BLS12_381_G1.mul(s % bls12_381.fields.Fr.ORDER),
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

const evalWithG1Setup2 = (p: Polynomial, setup: Setup2) => {
  const { coefficients } = p;
  return coefficients
    .map((coeff, index) => setup[index].g1.mul(coeff.n))
    .reduce((acc, curr) => acc.add(curr));
};

const p2 = new Polynomial(
  [1n, 2n, 3n, 4n, 7n, 7n, 7n, 7n, 13n, 13n, 13n, 13n, 13n, 13n, 13n, 13n],
  bls12_381.fields.Fr.ORDER
);

{
  const input = 1_927_409_816n;
  const e = p2.eval(input);
  const setup = trustedSetup(p2.degree() + 1, input);
  const setup2 = mytrustedSetup(p2.degree() + 1, input);

  for (let i = 0; i < setup.length; i++) {
    console.log(`==== s1: ${setup[i].s} ====`);
    console.log(setup[i].g1.x);
    console.log(`==== s2: ${setup2[i].s} ====`);
    console.log(setup2[i].g1.x.toString());
    console.log("================");
  }

  const e1 = evalWithG1Setup(p2, setup);
  const e2 = evalWithG1Setup2(p2, setup2);
  console.log(e1.x, e2.x.toString());
}
