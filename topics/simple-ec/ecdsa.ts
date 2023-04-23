import { print } from "../../utils/print";
import { EllipticCurve, Field, Point } from "./EllipticCurve";
import { random } from "./../random";

const BIG_NUM = 2n ** 256n;

const curve = EllipticCurve.SECP256K1;
const s = 0x83ecb3984a4f9ff03e84d5f9c0d7f888a81833643047acc58eb6431e01d9bac8n;
const S = curve.multiply(EllipticCurve.SECP256K1_G, s); //公開鍵
const msg = "こんにちは".split("").reduce((acc, c) => acc * 256n + BigInt(c.charCodeAt(0)), 0n);
const fr = new Field(0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141n); // 大きい素数

print("State", { s, S, msg });

const sign = (
  curve: EllipticCurve,
  G: Point,
  s: bigint,
  msg: bigint,
  k?: bigint
): [bigint, bigint] => {
  const z = fr.mod(msg);
  k ??= fr.mod(random(BIG_NUM)); //仮秘密鍵に相当する
  const Q = curve.multiply(G, k);
  const r = fr.mod(Q.x);
  return [r, fr.div(r * s + z, k)];
};

const signature = sign(curve, EllipticCurve.SECP256K1_G, s, msg);
print("Signature", signature);

const verify = (
  curve: EllipticCurve,
  G: Point,
  S: Point,
  msg: bigint,
  signature: [bigint, bigint]
) => {
  const z = fr.mod(msg);
  const [r, s] = signature;
  if (r === 0n || s === 0n) return false;

  const w = fr.inverse(s);
  const u1 = fr.mod(z * w); // == z / s
  const u2 = fr.mod(r * w); // == r / s
  const Q = curve.add(curve.multiply(G, u1), curve.multiply(S, u2));
  return r === fr.mod(Q.x);
};

print("Verify", { Ok: verify(curve, EllipticCurve.SECP256K1_G, S, msg, signature) });
