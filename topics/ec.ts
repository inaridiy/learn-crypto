import { EllipticCurve } from "./EllipticCurve";
import { print } from "../utils/print";

const curve = new EllipticCurve(4n, 14n);
const point1 = { x: -1n, y: 3n };
let current = point1;

for (let i = 1; i < 100; i++) {
  print(i.toString(), { ...current, OK: curve.isOnCurve(current) });
  current = curve.add(current, point1);
}
