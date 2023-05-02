import { ExtFQ } from "./ExtFQ";
import { FQ } from "./FQ";
import { Field } from "./types";

export type PointCoord<T extends Field<any, any> = Field<any, any>> = Readonly<{
  x: T;
  y: T;
}>;

export class EllipticCurve<T extends Field<unknown, unknown>> {
  constructor(a: T, b: T) {}
}
