import { EllipticCurve } from "./EllipticCurve";
import { print } from "../../utils/print";

const curve = EllipticCurve.SECP256K1;
const G = EllipticCurve.SECP256K1_G;

const q = curve.fp.prime;
const L = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141n;
const t = 1n + q - L;

// 1*P -t*P +q*P = O
const result = curve.add(
  curve.sub(curve.multiply(G, 1n), curve.multiply(G, t)),
  curve.multiply(G, q)
);

print("State1", { result });
