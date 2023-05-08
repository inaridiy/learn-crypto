import { pairing } from "./curve/pairing";
import { BLS12_381_G1, BLS12_381_G2, FQ12 } from "./curve/bls12-381";

const P1 = pairing(BLS12_381_G1, BLS12_381_G2);

console.log(P1.toString());
