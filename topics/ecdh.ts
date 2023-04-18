import { print } from "../utils/print";
import { EllipticCurve } from "./EllipticCurve";
import { random } from "./random";

const BIG_NUM = 2n ** 64n;

const curve = EllipticCurve.SECP256K1;
const a = random(BIG_NUM);
const b = random(BIG_NUM);

print("a", { a, b });

const aP = curve.multiply(EllipticCurve.SECP256K1_G, 7n);
// const bP = curve.multiply(EllipticCurve.SECP256K1_G, 12n);
