import mcl from "mcl-wasm";
import { print } from "../utils/print";

await mcl.init(mcl.BLS12_381);

const P = mcl.hashAndMapToG1("abc");
const Q = mcl.hashAndMapToG2("abc");

const n123 = new mcl.Fr();
n123.setInt(123);
const n456 = new mcl.Fr();
n456.setInt(456);

const aP = mcl.mul(P, n123);
const bQ = mcl.mul(Q, n456);

const e1 = mcl.pairing(P, Q); // e(P, Q)
const e2 = mcl.pairing(aP, bQ); // e(aP, bQ)

const e3 = mcl.pow(e1, mcl.mul(n123, n456)); // e(P, Q)^(n123 * n456)

print("Pairing", {
  e1,
  e2,
  e3,
  Ok: e2.isEqual(e3), // e(aP, bQ) == e(P, Q)^(n123 * n456)
});
