export type FieldFactory<T, TLike> = {
  zero(): T;
  one(): T;
  from(value: bigint): T;
  from(value: TLike): T;
  from(value: T): T;
};

export interface Field<T, TLike> {
  one(): Field<T, TLike>;
  zero(): Field<T, TLike>;
  clone(): Field<T, TLike>;
  extend(value: this | TLike): Field<T, TLike>;
  eq(other: this | TLike): boolean;
  isZero(): boolean;

  add(other: this | TLike): Field<T, TLike>;
  sub(other: this | TLike): Field<T, TLike>;
  mul(other: this | TLike): Field<T, TLike>;

  div(other: this | TLike): Field<T, TLike>;

  pow(n: bigint): Field<T, TLike>;
}
