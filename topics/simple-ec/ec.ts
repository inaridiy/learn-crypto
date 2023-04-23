import { EllipticCurve } from "./EllipticCurve";
import { print } from "../../utils/print";

const curve = EllipticCurve.SECP256K1;
const point1 = EllipticCurve.SECP256K1_G;
let current = EllipticCurve.ZERO_POINT;

console.time("add");

for (let i = 0; i < 1000; i++) {
  current = curve.add(current, point1);
}

const ans = curve.multiply(point1, 1000n);
print("Result", { current, ans });
console.timeEnd("add");
