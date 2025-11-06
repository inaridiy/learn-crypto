import { CommutativeRingElement, EuclideanRingElement } from "./interface";

export const modBigint = (a: bigint, m: bigint) => {
  const res = a % m;
  return res >= 0n ? res : res + m;
};

export const pow = <T extends CommutativeRingElement<T, TLike>, TLike = T>(field: T, n: bigint) => {
  let result = field.structure.one();
  let base = field.clone();
  while (n > 0n) {
    if (n % 2n === 1n) result = result.mul(base);
    if (n > 0n) base = base.mul(base);
    n >>= 1n;
  }
  return result;
};

export const extendedGCD = <T extends EuclideanRingElement<T, TLike>, TLike = T>(a: T, b: T) => {
  let [r1, r2] = [a, b];
  let [s1, s2] = [a.structure.one(), a.structure.zero()];
  let [t1, t2] = [a.structure.zero(), a.structure.one()];

  while (!r2.isZero()) {
    const [q, r] = r1.divmod(r2);
    [r1, r2] = [r2, r];
    [s1, s2] = [s2, s1.sub(q.mul(s2))];
    [t1, t2] = [t2, t1.sub(q.mul(t2))];
  }

  return [r1, s1, t1];
};

export const extendedGCDBigint = (a: bigint, b: bigint) => {
  let [r1, r2] = [a, b];
  let [s1, s2] = [1n, 0n];
  let [t1, t2] = [0n, 1n];

  while (r2 !== 0n) {
    const q = r1 / r2;
    [r1, r2] = [r2, r1 - q * r2];
    [s1, s2] = [s2, s1 - q * s2];
    [t1, t2] = [t2, t1 - q * t2];
  }

  return [r1, s1, t1];
};
