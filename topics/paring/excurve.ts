import { ExtFQ, ExtFQFactory, ExtFQLike } from "./curve/ExtFQ";
import { CurvePoint, EllipticCurve } from "./curve/EllipticCurve";

const q =
  4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787n;

const FQ = new ExtFQFactory(q, 1n);
const FQCurve = new EllipticCurve(0n, 4n, FQ);
const G1 = new CurvePoint(
  3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507n,
  1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569n,
  FQCurve
);

const FQ2 = new ExtFQFactory(q, [1n, 0n, 1n]);
const FQ2Curve = new EllipticCurve(FQ2.from(0n), FQ2.from([4n, 4n]), FQ2);
const G2 = new CurvePoint(
  FQ2.from([
    352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160n,
    3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758n,
  ]),
  FQ2.from([
    1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905n,
    927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582n,
  ]),
  FQ2Curve
);

const P3 = G1.mul(100n);

console.log(P3.x.toString());