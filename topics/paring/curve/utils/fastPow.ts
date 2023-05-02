import { Field } from "../types";

export const fastPow = <T extends Field<any, any>>(filed: T, n: bigint) => {
  if (n === 0n) return filed.zero();

  let result: T = filed.one();
  let base = filed.clone();
  while (n > 0n) {
    if (n % 2n === 1n) result = result.mul(base);
    if (n > 0n) base = base.mul(base);
    n >>= 1n;
  }

  return result;
};
