import { print } from "../../../../utils/print";
import { Polynomial } from "../Polynomial";

export const extGCD = <T = bigint>(a: T, b: T): [T, T, T] => {
  if (typeof a === "bigint") return bigintExtGCD(a, b as bigint) as any;
  else if (a instanceof Polynomial) return polyExtGCD(a, b as Polynomial) as any;
  else throw new Error("Not implemented");
};

export const bigintExtGCD = (a: bigint, b: bigint): [bigint, bigint, bigint] => {
  let x = 0n;
  let y = 1n;
  let u = 1n;
  let v = 0n;
  while (a !== 0n) {
    const q = b / a;
    const r = b % a;
    const m = x - u * q;
    const n = y - v * q;
    b = a;
    a = r;
    x = u;
    y = v;
    u = m;
    v = n;
  }
  return [b, x, y];
};

//上と実装が違うが、試行錯誤していたため
export const polyExtGCD = (a: Polynomial, b: Polynomial): [Polynomial, Polynomial, Polynomial] => {
  if (a.p !== b.p) throw new Error("Must be same field");

  let [r1, r2] = [a, b];
  let [s1, s2] = [new Polynomial([1n], a.p), new Polynomial([0n], a.p)];
  let [t1, t2] = [new Polynomial([0n], a.p), new Polynomial([1n], a.p)];

  while (!r2.isZero()) {
    const q = r1.div(r2);
    [r1, r2] = [r2, r1.sub(q.mul(r2))];
    [s1, s2] = [s2, s1.sub(q.mul(s2))];
    [t1, t2] = [t2, t1.sub(q.mul(t2))];
  }

  return [r1, s1, t1];
};
