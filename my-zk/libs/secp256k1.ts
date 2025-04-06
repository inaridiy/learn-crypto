import { EllipticCurve } from "./curve";
import { FiniteFieldFactory } from "./finite-field";

export const FQ = new FiniteFieldFactory(2n ** 256n - 2n ** 32n - 977n);

export const SECP256K1_G = {
  x: FQ.from(55066263022277343669578718895168534326250603453777594175500187360389116729240n),
  y: FQ.from(32670510020758816978083085130507043184471273380659243275938904335757337482424n),
};

export const SECP256K1 = new EllipticCurve(FQ.from(0n), FQ.from(7n));
