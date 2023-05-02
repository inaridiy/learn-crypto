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

export type FieldFactory<T, TLike> = {
  zero(): Field<T, TLike>;
  one(): Field<T, TLike>;
  from(value: bigint): Field<T, TLike>;
  from(value: TLike): Field<T, TLike>;
};
