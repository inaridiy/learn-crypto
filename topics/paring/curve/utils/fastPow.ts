import { Field } from "../types";

export const fastPow = <T>(filed: Field<T, any>, n: bigint): T => {
  if (n === 0n) return filed.zero() as T;

  let result = filed.one();
  let base = filed.clone();
  while (n > 0n) {
    if (n % 2n === 1n) result = result.mul(base);
    if (n > 0n) base = base.mul(base);
    n >>= 1n;
  }

  return result as T;
};
