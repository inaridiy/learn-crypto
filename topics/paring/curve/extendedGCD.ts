import { Polynomial } from "./Polynomial";

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

export const polyExtGCD = (a: Polynomial, b: Polynomial): [Polynomial, Polynomial, Polynomial] => {
  if (a.p !== b.p) throw new Error("Must be same field");
  let x = Polynomial.zero(a.p);
  let y = new Polynomial([1n], a.p);
  let u = new Polynomial([1n], a.p);
  let v = Polynomial.zero(a.p);

  while (!a.isZero()) {
    const q = b.div(a);
    const r = b.mod(a);
    const m = x.sub(u.mul(q));
    const n = y.sub(v.mul(q));
    b = a;
    a = r;
    x = u;
    y = v;
    u = m;
    v = n;
  }

  return [b, x, y];
};
