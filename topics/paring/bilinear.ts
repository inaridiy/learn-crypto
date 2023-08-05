import { pairing } from "./curve/pairing";
import { BLS12_381_G1, BLS12_381_G2 } from "./curve/bls12-381";

if (import.meta.vitest) {
  it("pairing", async () => {
    const a = 123n;
    const b = 456n;

    const aP = BLS12_381_G1.mul(a);
    const bQ = BLS12_381_G2.mul(b);

    const e1 = pairing(BLS12_381_G1, BLS12_381_G2);
    const e2 = pairing(aP, bQ);

    expect(e1.pow(a * b).value.toString()).toEqual(e2.value.toString()); //双線形だぜ☆（＾～＾）
  });
}
