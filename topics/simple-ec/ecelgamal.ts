import { EllipticCurve, Point } from "./EllipticCurve";
import { print } from "../../utils/print";
import { random } from "../random";
import { sqrtMod } from "../math";

const curve = EllipticCurve.SECP256K1;
const G = EllipticCurve.SECP256K1_G;

const encode = (msg: bigint) => {
  const L = 16n;
  for (let i = 0n; i < 2n ** L; i++) {
    const x = (msg << L) | i;
    const y2 = curve.calcY2(x);
    const [isSquare, y] = sqrtMod(y2, EllipticCurve.BIG_PRIME);
    if (isSquare) return { x, y };
  }
  throw new Error("Failed to encode");
};

const decode = (p: Point) => p.x >> 16n;

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

const msg = 0x123456789n;
const encoded = encode(msg);

print("State", { s, sP, msg });

const [c1, c2] = encrypt(sP, encoded);

print("Encrypted", { c1, c2 });

const dec = decrypt(s, c1, c2);
const decMsg = decode(dec);

print("Decrypted", { dec, decMsg });
