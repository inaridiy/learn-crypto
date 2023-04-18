import { print } from "../utils/print";
import { EllipticCurve } from "./EllipticCurve";
import { random } from "./random";

const BIG_NUM = 2n ** 256n;

const curve = EllipticCurve.SECP256K1;
const a = random(BIG_NUM);
const b = random(BIG_NUM);

print("Private", { a, b });

const aP = curve.multiply(EllipticCurve.SECP256K1_G, a);
const bP = curve.multiply(EllipticCurve.SECP256K1_G, b);

print("Public", { aP, bP });

const baP = curve.multiply(aP, b);
const abP = curve.multiply(bP, a);

print("Shared", { baP, abP, OK: baP.x === abP.x });
