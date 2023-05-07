import { ExtFQFactory } from "./ExtFQ";

const q =
  4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787n;

export const FQ = new ExtFQFactory(q);
export const FQ2 = new ExtFQFactory(q, [1n, 0n, 1n]);
export const FQ12 = new ExtFQFactory(q, [2n, 0n, 0n, 0n, 0n, 0n, -2n, 0n, 0n, 0n, 0n, 0n, 1n]);
