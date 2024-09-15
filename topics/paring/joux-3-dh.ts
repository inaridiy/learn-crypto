import { CurvePoint } from "./curve/EllipticCurve";
import { pairing } from "./curve/pairing";
import { BLS12_381_G1, BLS12_381_G2 } from "./curve/bls12-381";
import { print } from "../../utils/print";

const joux = (sk: bigint, bG1: CurvePoint, cG2: CurvePoint) => {
  return pairing(bG1, cG2).pow(sk);
};

if (import.meta.vitest) {
  it("joux-3-dh", async () => {
    const a = 1234n;
    const b = 5678n;
    const c = 9012n;

    const aG1 = BLS12_381_G1.mul(a);
    const aG2 = BLS12_381_G2.mul(a);

    const bG1 = BLS12_381_G1.mul(b);
    const bG2 = BLS12_381_G2.mul(b);

    const cG1 = BLS12_381_G1.mul(c);
    const cG2 = BLS12_381_G2.mul(c);

    const sharedA = joux(a, bG1, cG2);
    const sharedB = joux(b, cG1, aG2);
    const sharedC = joux(c, aG1, bG2);

    print("shared Keys", {
      a: sharedA.value.toString(),
      b: sharedB.value.toString(),
      c: sharedC.value.toString(),
    });
    expect(sharedA.value.toString()).toEqual(sharedB.value.toString());
    expect(sharedA.value.toString()).toEqual(sharedC.value.toString());
  });
}
