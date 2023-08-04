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

if (import.meta.vitest) {
  it("lifted ec elgamal", async () => {
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

    expect(decMsg).toEqual(msg);
  });

  it("additive quasi-homomorphism of lifted ec elgamal", async () => {
    const s = random(EllipticCurve.BIG_PRIME - 1n);
    const sP = curve.multiply(G, s);

    const msg1 = 0x123n;
    const encoded1 = encode(msg1);

    const msg2 = 0x456n;
    const encoded2 = encode(msg2);

    print("State", { s, sP, msg1, msg2 });

    const encrypted1 = encrypt(sP, encoded1);
    const encrypted2 = encrypt(sP, encoded2);

    print("Encrypted", { encrypted1, encrypted2 });

    const c1Sum = curve.add(encrypted1[0], encrypted2[0]);
    const c2Sum = curve.add(encrypted1[1], encrypted2[1]);

    print("Sum", { c1Sum, c2Sum });

    const dec = decrypt(s, c1Sum, c2Sum);
    const decMsg = decode(dec);

    print("Decrypted", { dec, decMsg, Ok: decMsg === msg1 + msg2 });

    expect(decMsg).toEqual(msg1 + msg2);
  });
}
