import { EllipticCurve, Point } from "./EllipticCurve";
import { print } from "../../utils/print";
import { random } from "../random";

const curve = EllipticCurve.SECP256K1;
const G = EllipticCurve.SECP256K1_G;

const encode = (msg: bigint) => {
  return curve.multiply(G, msg);
};

const decode = (p: Point) => {
  const MaxMsgSize = 2n ** 32n;
  for (let i = 0n; i < MaxMsgSize; i++) {
    if (curve.multiply(G, i).x === p.x) return i;
  }
  throw new Error("Failed to decode");
};

const encrypt = (sP: Point, m: Point) => {
  const t = random(EllipticCurve.BIG_PRIME - 1n);
  const c1 = curve.add(m, curve.multiply(sP, t));
  const c2 = curve.multiply(G, t);
  return [c1, c2] as const;
};

const decrypt = (s: bigint, c1: Point, c2: Point) => {
  const dec = curve.sub(c1, curve.multiply(c2, s));
  return dec;
};

const s = random(EllipticCurve.BIG_PRIME - 1n);
const sP = curve.multiply(G, s);

const msg = 0x1234n;
const encoded = encode(msg);

print("State", { s, sP, msg });

const [c1, c2] = encrypt(sP, encoded);

print("Encrypted", { c1, c2 });

const dec = decrypt(s, c1, c2);
const decMsg = decode(dec);

print("Decrypted", { dec, decMsg });
