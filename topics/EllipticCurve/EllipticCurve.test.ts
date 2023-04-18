import { EllipticCurve } from "./EllipticCurve";

it("Test Curve", () => {
  const curve = new EllipticCurve(4n, 14n);
  const point1 = { x: -1n, y: 3n };
  expect(curve.isOnCurve(point1)).toBe(true);
  expect(curve.isOnCurve({ x: 0n, y: 0n })).toBe(false);
  expect(curve.isOnCurve(curve.add(point1, point1))).toBe(true);
  expect(curve.isOnCurve(curve.add(point1, point1, point1, point1))).toBe(true);
  // expect(curve.multiply(point1, 16n)).toEqual(curve.add(...Array(16).fill(point1)));
});

it("Test SECP256K1", () => {
  const curve = EllipticCurve.SECP256K1;
  const point = EllipticCurve.SECP256K1_G;

  expect(curve.isOnCurve(point)).toBe(true);
});
