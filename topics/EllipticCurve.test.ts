import { EllipticCurve } from "./EllipticCurve";

it("Test Curve", () => {
  const curve = new EllipticCurve(4n, 14n);
  const point1 = { x: -1n, y: 3n };
  expect(curve.isOnCurve(point1)).toBe(true);

  expect(curve.isOnCurve(curve.add(point1, point1))).toBe(true);
  expect(curve.isOnCurve(curve.add(point1, point1, point1))).toBe(true);
  console.log(curve.multiply(point1, 3n));
  expect(curve.add(point1, point1, point1)).toEqual(curve.multiply(point1, 3n));
});
