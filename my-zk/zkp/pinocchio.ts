import { random } from "../../topics/random";
import { BLS12_381_G1, BLS12_381_G2 } from "../curves/bls12-381";
import { FiniteFieldElement } from "../primitive/finite-field";
import { matrixToPolynomials } from "../primitive/lagrange";
import { R1CSConstraints } from "./r1cs";

export const trustedSetup = (r1cs: R1CSConstraints<FiniteFieldElement>) => {
  const [aps, bps, cps] = [r1cs.a, r1cs.b, r1cs.c].map((v) =>
    matrixToPolynomials(r1cs.structure, v)
  );
  const maxDgree = Math.max(...[...aps, ...bps, ...cps].map((p) => p.degree()));

  const [s, alpha, betaA, betaB, betaC, gamma] = Array.from({ length: 6 }, () =>
    random(r1cs.structure.p)
  ).map((v) => r1cs.structure.from(v));

  const [g1, g2] = [BLS12_381_G1.generator(), BLS12_381_G2.generator()];

  const ss = Array.from({ length: maxDgree }, (_, i) => s.pow(BigInt(i + 1)));

  const setup = (v: FiniteFieldElement) => ({ g1: g1.scale(v.n), g2: g2.scale(v.n) });

  const tOne = setup(r1cs.structure.one());
  const tAlpha = setup(alpha);
  const tGamma = setup(gamma);
  const tBetaAgamma = setup(gamma.mul(betaA));
  const tBetaBgamma = setup(gamma.mul(betaB));
  const tBetaCgamma = setup(gamma.mul(betaC));

  const ts = ss.map((v) => setup(v));
  const tAlphaS = ss.map((v) => alpha.mul(v)).map((v) => setup(v));

  const tas = aps.map((p) => p.eval(s)).map((v) => setup(v));
  const tbs = bps.map((p) => p.eval(s)).map((v) => setup(v));
  const tcs = cps.map((p) => p.eval(s)).map((v) => setup(v));

  const tBetaAas = aps.map((p) => p.eval(s).mul(betaA)).map((v) => setup(v));
  const tBetaBbs = bps.map((p) => p.eval(s).mul(betaB)).map((v) => setup(v));
  const tBetaCcs = cps.map((p) => p.eval(s).mul(betaC)).map((v) => setup(v));

  const midStart = r1cs.index.output + 1;

  const ek = {
    ts,
    tAlphaS,
    tasMid: tas.slice(midStart),
    tbsMid: tas.slice(midStart),
    tcsMid: tas.slice(midStart),
    tBetaAasMid: tBetaAas.slice(midStart),
    tBetaBbsMid: tBetaBbs.slice(midStart),
    tBetaCcsMid: tBetaCcs.slice(midStart),
  };

  const vk = {
    tOne,
    tAlpha,
    tGamma,
    tBetaAgamma,
    tBetaBgamma,
    tBetaCgamma,
    //...
  };
};
