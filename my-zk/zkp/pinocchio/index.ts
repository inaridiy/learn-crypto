import { random } from "../../../topics/random";
import { BLS12_381_G1, BLS12_381_G2 } from "../../primitive/curves/bls12-381";
import { pairing } from "../../primitive/paring";
import {
  createProve,
  createTrustedSetup,
  createVerify,
  createZkProve,
  createZkTrustedSetup,
} from "./protocol";
import { createPointOps } from "./operations";
import { PinocchioConfig, PinocchioInstance, PinocchioRuntime } from "./types";

export * from "./types";

const createRuntime = (config: PinocchioConfig): PinocchioRuntime => ({
  pairing: config.pairing,
  sampleFieldElement: config.sampleFieldElement,
  pointOps: createPointOps(config.generators),
});

/**
 * pinocchio関数:
 * 外部から有限体・ペアリング群などの数学的構造を注入し、
 * その構造上で動作する trustedSetup / prove / verify のセットを生成する。
 */
export const pinocchio = (config: PinocchioConfig): PinocchioInstance => {
  const runtime = createRuntime(config);
  const trustedSetup = createTrustedSetup(runtime);
  const prove = createProve(runtime);

  const zkTrustedSetup = createZkTrustedSetup(runtime);
  const zkProve = createZkProve(runtime);

  const verify = createVerify(runtime);

  return { trustedSetup, prove, zkTrustedSetup, zkProve, verify };
};

export const defaultPinocchioConfig: PinocchioConfig = {
  generators: {
    g1: BLS12_381_G1.generator(),
    g2: BLS12_381_G2.generator(),
  },
  pairing,
  sampleFieldElement: (field) => field.from(random(field.p)),
};
