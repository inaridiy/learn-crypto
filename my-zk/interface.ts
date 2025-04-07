export interface FieldFactory<T extends Field<T, FieldLike>, FieldLike> {
  zero(): T;
  one(): T;
  from(value: FieldLike): T;
}

export interface Field<T extends Field<T, FieldLike>, FieldLike> {
  factory: FieldFactory<T, FieldLike>;
  eq(other: T): boolean;
  isZero(): boolean;
  clone(): T;
  add(other: T): T;
  sub(other: T): T;
  mul(other: T): T;
  scale(n: bigint): T;
  div(other: T): T;
  pow(n: bigint): T;
}
