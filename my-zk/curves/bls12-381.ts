import { ECPointCyclicGroup, EllipticCurve } from "../primitive/elliptic-curve";
import { FiniteField } from "../primitive/finite-field";
import { PolynomialFactory, type PolynomialOnFF } from "../primitive/polynomial";
import {
  ExtendedFiniteField,
  ExtendedFiniteFieldElement,
} from "../primitive/extended-finite-field";

export const FIELD_MODULUS =
  4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787n;

export const FQ = new FiniteField(FIELD_MODULUS);
export const POLY = new PolynomialFactory(FQ);

export const CURVE_ORDER =
  52435875175126190479447740508185965837690552500527637822603658699938581184513n;

export const FQ1 = new ExtendedFiniteField(POLY.from([0n, 1n]));

export const FQ2 = new ExtendedFiniteField(POLY.from([1n, 0n, 1n]));

export const FQ12 = new ExtendedFiniteField(
  POLY.from([2n, 0n, 0n, 0n, 0n, 0n, -2n, 0n, 0n, 0n, 0n, 0n, 1n])
);

export const BLS12_381_G1 = new ECPointCyclicGroup<ExtendedFiniteFieldElement, PolynomialOnFF>(
  FQ1.zero(),
  FQ1.from(POLY.from([4n])),
  {
    x: FQ1.from(
      POLY.from([
        3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507n,
      ])
    ),
    y: FQ1.from(
      POLY.from([
        1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569n,
      ])
    ),
  },
  10n
);
export const BLS12_381_G2 = new ECPointCyclicGroup<ExtendedFiniteFieldElement, PolynomialOnFF>(
  FQ2.zero(),
  FQ2.from(POLY.from([4n, 4n])),
  {
    x: FQ2.from(
      POLY.from([
        352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160n,
        3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758n,
      ])
    ),
    y: FQ2.from(
      POLY.from([
        1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905n,
        927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582n,
      ])
    ),
  },
  10n
);
export const BLS12_381_FQ12 = new EllipticCurve<ExtendedFiniteFieldElement, PolynomialOnFF>(
  FQ12.zero(),
  FQ12.from(POLY.from([4n]))
);
