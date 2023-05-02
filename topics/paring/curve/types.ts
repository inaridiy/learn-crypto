export interface Field<T, TLike> {
  one(): T;
  zero(): T;
  clone(): T;
  extend(value: TLike): T;
  eq(other: TLike): boolean;
  isZero(): boolean;

  add(other: TLike): T;
  sub(other: TLike): T;
  mul(other: TLike): T;

  div(other: TLike): T;

  pow(n: bigint): T;
}

export type FieldFactory<T, TLike> = {
  zero(): T;
  one(): T;
  from(value: TLike): T;
};
