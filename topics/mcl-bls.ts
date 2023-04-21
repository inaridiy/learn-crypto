import mcl from "mcl-wasm";
import { print } from "../utils/print";

await mcl.init(mcl.BLS12_381);

const Q = mcl.hashAndMapToG2("abc");

const sign = (msg: string, sk: mcl.Fr): mcl.G1 => {
  const hash = mcl.hashAndMapToG1(msg);
  return mcl.mul(hash, sk);
};

const verify = (msg: string, sig: mcl.G1, pk: mcl.G2): boolean => {
  const hash = mcl.hashAndMapToG1(msg);
  const e1 = mcl.pairing(sig, Q);
  const e2 = mcl.pairing(hash, pk);

  return e1.isEqual(e2);
};

const msg = "Hello World";
const sk = new mcl.Fr();
sk.setByCSPRNG();
const pk = mcl.mul(Q, sk);

print("State", { msg, sk: sk.serializeToHexStr(), pk: pk.serializeToHexStr() });

const sig = sign(msg, sk);
const valid = verify(msg, sig, pk);

print("Valid", { sig: sig.serializeToHexStr(), Ok: valid });
