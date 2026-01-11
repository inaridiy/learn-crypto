import { BLS12_381_G2, BLS12_381_G1 } from "./primitive/curves/bls12-381";
import { ECPointOnExtendedFF } from "./primitive/elliptic-curve";
import { pairing } from "./primitive/paring";

const s = 0x83ecb3984a4f9ff03e84d5f9c0d7f888a81833643047acc58eb6431e01d9bac8n;
const sG = BLS12_381_G2.generator().scale(s); //公開鍵

const hash = (msg: string) => {
  const hashedBigInt = BigInt(
    "0x" + msg.split("").reduce((acc, c) => acc + c.charCodeAt(0).toString(16), "")
  );
  return BLS12_381_G1.generator().scale(hashedBigInt);
};

export const sign = (msg: string, s: bigint) => {
  const hashed = hash(msg);
  const sig = hashed.scale(s);
  return sig;
};

export const verify = (msg: string, sG: ECPointOnExtendedFF, sig: ECPointOnExtendedFF) => {
  const hashed = hash(msg);
  const e1 = pairing(sig, BLS12_381_G2.generator()); //e(hash(msg)*s, G2)
  const e2 = pairing(hashed, sG); //e(hash(msg), s*G2)
  return e1.eq(e2);
};

if (import.meta.vitest) {
  it("BLS", async () => {
    const msg = "hello world";
    const sig = sign(msg, s);
    const verified = verify(msg, sG, sig);

    expect(verified).toEqual(true);
  });
}
