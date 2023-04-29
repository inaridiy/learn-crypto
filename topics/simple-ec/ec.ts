import { EllipticCurve } from "./EllipticCurve";
import { print } from "../../utils/print";

const curve = EllipticCurve.SECP256K1;
const point1 = EllipticCurve.SECP256K1_G;
const zero = EllipticCurve.ZERO_POINT;

const L = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141n;

const zero2 = curve.multiply(point1, L + 1n); // 0n
print(zero2);
