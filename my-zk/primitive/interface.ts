// 可換環（Commutative Ring）
export interface CommutativeRing<
  TElement extends CommutativeRingElement<TElement, TLike>,
  TLike = unknown
> {
  zero: () => TElement;
  one: () => TElement;

  from: (value: TElement | TLike) => TElement;
}

// 可換環の元（Commutative　Ring element）
export interface CommutativeRingElement<
  TSelf extends CommutativeRingElement<TSelf, TLike>,
  TLike = unknown
> {
  readonly structure: CommutativeRing<TSelf, TLike>;
  readonly CommutativeRingElement: true;

  clone: () => TSelf;

  eq: (other: TSelf) => boolean;
  isZero: () => boolean;
  isOne: () => boolean;

  add: (other: TSelf) => TSelf;
  sub: (other: TSelf) => TSelf;
  mul: (other: TSelf) => TSelf;

  pow: (n: bigint) => TSelf;
  scale: (n: bigint) => TSelf;
}

export interface Field<TElement extends FieldElement<TElement, TLike>, TLike = unknown>
  extends CommutativeRing<TElement, TLike> {}

export interface FieldElement<TSelf extends FieldElement<TSelf, TLike>, TLike = unknown>
  extends CommutativeRingElement<TSelf, TLike> {
  readonly FieldElement: true;
  inverse: () => TSelf;
  div: (other: TSelf) => TSelf;
}

//ガバガバだけど、ユークリッドの互除法が使える環
export interface EuclidRingElement<
  TSelf extends CommutativeRingElement<TSelf, TLike>,
  TLike = unknown
> extends CommutativeRingElement<TSelf, TLike> {
  readonly EuclidRingElement: true;
  div: (other: TSelf) => TSelf;
}
