import { EllipticCurve } from "./EllipticCurve";
import { print } from "../utils/print";

const curve = new EllipticCurve(4n, 14n);
const point1 = { x: -1n, y: 3n };
let current = point1;

for (let i = 1; i < 16; i++) {
  print(i.toString(), { ...current, OK: curve.isOnCurve(current) });
  current = curve.add(current, point1);
}

print("point1", { ...curve.multiply(point1, 16n), OK: curve.isOnCurve(point1) });