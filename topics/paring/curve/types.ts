export type FieldFactory<T, TLike> = {
  zero(): T;
  one(): T;
  from(value: bigint): T;
  from(value: TLike): T;
  from(value: T): T;
};
