export const fastPow = (x: bigint, n: bigint, mod: bigint) => {
  let res = 1n;
  while (n > 0n) {
    if (n & 1n) res = (res * x) % mod;
    x = (x * x) % mod;
    n >>= 1n;
  }
  return res;
};

export const modInverse = (a: bigint, mod: bigint) => {
  return fastPow(a, mod - 2n, mod);
};

export const sqrtMod = (a: bigint, p: bigint): [boolean, bigint] => {
  if (p % 4n !== 3n) throw new Error("p%4 must be 3");
  if (fastPow(a, (p - 1n) / 2n, p) !== 1n) return [false, 0n]; // no solution
  const root = fastPow(a, (p + 1n) / 4n, p);
  return [true, root];
};
