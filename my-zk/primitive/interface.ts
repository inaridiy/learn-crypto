// =================================================================
// 1. 基本的なインターフェース
// =================================================================

/**
 * 何らかの代数構造を持つ元の基本インターフェース
 * @template TSelf - 自分自身の型 (CRTP: Curiously Recurring Template Pattern)
 * @template TLike - 自分自身に変換可能な、より単純な型 (例: number, bigint)
 */
export interface Element<TSelf extends Element<TSelf, TLike>, TLike = unknown> {
  /** この元が属する代数構造への参照 */
  readonly structure: AlgebraicStructure<TSelf, TLike>;

  /** 新しいインスタンスとして自身を複製する */
  clone: () => TSelf;

  /** 他の元と値が等しいか判定する */
  eq: (other: TSelf) => boolean;
}

/**
 * 代数構造そのものを表す基本インターフェース
 */
export interface AlgebraicStructure<TElement extends Element<TElement, TLike>, TLike = unknown> {
  /** プリミティブな値や他のインスタンスから、この構造の元を生成する */
  from: (value: TElement | TLike) => TElement;
}

// =================================================================
// 2. 群の階層 (環の加法の基礎となる部分)
// =================================================================

/**
 * アーベル群（可換群）を表す構造。環の加法群のモデルとなる。
 */
export interface AbelianGroup<
  TElement extends AbelianGroupElement<TElement, TLike>,
  TLike = unknown
> extends AlgebraicStructure<TElement, TLike> {
  /** 加法の単位元 (ゼロ元) */
  zero: () => TElement;
}

/**
 * アーベル群の元
 */
export interface AbelianGroupElement<
  TSelf extends AbelianGroupElement<TSelf, TLike>,
  TLike = unknown
> extends Element<TSelf, TLike> {
  readonly structure: AbelianGroup<TSelf, TLike>;

  isZero: () => boolean;
  add: (other: TSelf) => TSelf;
  sub: (other: TSelf) => TSelf;
  /** 加法の逆元（反数）: -this */
  negate: () => TSelf;
  /** スカラー倍（整数倍）: n * this */
  scale: (n: bigint) => TSelf;
}

/**
 * 巡回群を表す構造
 */
export interface CyclicGroup<TElement extends AbelianGroupElement<TElement, TLike>, TLike = unknown>
  extends AbelianGroup<TElement, TLike> {
  /** 群の位数 (order)。無限巡回群の場合は "inf" */
  readonly order: bigint | "inf";
  /** 群の生成元 */
  generator: () => TElement;
}

// =================================================================
// 3. 環の階層
// =================================================================

/**
 * (単位元を持つ)可換環を表す構造
 * 加法についてはアーベル群の性質を継承する
 */
export interface CommutativeRing<
  TElement extends CommutativeRingElement<TElement, TLike>,
  TLike = unknown
> extends AbelianGroup<TElement, TLike> {
  /** 乗法の単位元 (1) */
  one: () => TElement;
}

/**
 * (単位元を持つ)可換環の元
 */
export interface CommutativeRingElement<
  TSelf extends CommutativeRingElement<TSelf, TLike>,
  TLike = unknown
> extends AbelianGroupElement<TSelf, TLike> {
  readonly structure: CommutativeRing<TSelf, TLike>;

  isOne: () => boolean;
  mul: (other: TSelf) => TSelf;
  /** べき乗: this ^ n (nは非負整数) */
  pow: (n: bigint) => TSelf;
}

/**
 * ユークリッド環を表す構造
 */
export interface EuclideanRing<
  TElement extends EuclideanRingElement<TElement, TLike>,
  TLike = unknown
> extends CommutativeRing<TElement, TLike> {}

/**
 * ユークリッド環の元。ユークリッド除算（余りのある割り算）が可能。
 */
export interface EuclideanRingElement<
  TSelf extends EuclideanRingElement<TSelf, TLike>,
  TLike = unknown
> extends CommutativeRingElement<TSelf, TLike> {
  readonly structure: EuclideanRing<TSelf, TLike>;

  /** 商 (Quotient): this / other の整数商 */
  quotient: (other: TSelf) => TSelf;

  /** 剰余 (Remainder): this % other */
  remainder: (other: TSelf) => TSelf;

  /** 商と剰余を同時に計算する */
  divmod: (other: TSelf) => readonly [quotient: TSelf, remainder: TSelf];
}

/**
 * 体 (Field) を表す構造。体はユークリッド環でもある。
 */
export interface Field<TElement extends FieldElement<TElement, TLike>, TLike = unknown>
  extends EuclideanRing<TElement, TLike> {}

/**
 * 体の元。ゼロ元以外で除算が可能。
 */
export interface FieldElement<TSelf extends FieldElement<TSelf, TLike>, TLike = unknown>
  extends EuclideanRingElement<TSelf, TLike> {
  readonly structure: Field<TSelf, TLike>;

  /** 乗法の逆元: 1 / this */
  inverse: () => TSelf;

  /** 除算: this / other */
  div: (other: TSelf) => TSelf;
}
