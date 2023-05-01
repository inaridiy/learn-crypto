import { FQ } from "./curve/FQ";
import { ExtFQ } from "./curve/ExtFQ";
import { Polynomial } from "./curve/Polynomial";
import { print } from "../../utils/print";

const q = 2n;
const extDegree = 2;
const extModPoly = [1n, 1n, 1n];

const result1: bigint[][] = [];
for (let i = 0n; i < q; i++) {
  result1[Number(i)] = [];
  for (let j = 0n; j < q; j++) {
    const a = new FQ(q, i);
    const b = new FQ(q, j);
    const c = a.mul(b);
    result1[Number(i)][Number(j)] = c.n;
  }
}

console.table(result1);

const extFQ = new ExtFQ(q, extModPoly);
const values = [[0n], [1n], [0n, 1n], [1n, 1n]]; // 0, 1, x, x + 1
const result2: string[][] = [];
for (const i in values) {
  result2[Number(i)] = [];
  for (const j in values) {
    const a = extFQ.mul(values[i], values[j]);
    result2[Number(i)][Number(j)] = a.toString();
  }
}

console.table(result2);
