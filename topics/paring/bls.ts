import { CurvePoint } from "./curve/EllipticCurve";
import { pairing } from "./curve/pairing";
import { BLS12_381_G1, BLS12_381_G2 } from "./curve/bls12-381";
import { print } from "../../utils/print";

const s = 0x83ecb3984a4f9ff03e84d5f9c0d7f888a81833643047acc58eb6431e01d9bac8n;
const Q = BLS12_381_G2.mul(s); //公開鍵

const hash = (msg: string) => {
  const hashedBigInt = BigInt(
    "0x" + msg.split("").reduce((acc, c) => acc + c.charCodeAt(0).toString(16), "")
  );
  return BLS12_381_G1.mul(hashedBigInt);
};

const sign = (msg: string) => {
  const hashed = hash(msg);
  const sig = hashed.mul(s);
  return sig;
};

const verify = (msg: string, sig: CurvePoint) => {
  const hashed = hash(msg);
  const e1 = pairing(sig, BLS12_381_G2); //e(hash(msg)*s, G2)
  const e2 = pairing(hashed, Q); //e(hash(msg), s*G2)
  return e1.eq(e2);
};

const msg = "hello world";
const sig = sign(msg);
print("Verify", { Ok: verify(msg, sig) });
