import { Field, FieldFactory } from "./interface";
import { Polynomial } from "./polynomial";

export class ExtendedFiniteFieldFactory implements FieldFactory<ExtendedFiniteField, Polynomial> {
  constructor(public readonly mod: Polynomial) {}

  zero(): ExtendedFiniteField {
    return new ExtendedFiniteField(this, this.mod, this.mod.factory.zero());
  }

  one(): ExtendedFiniteField {
    return new ExtendedFiniteField(this, this.mod, this.mod.factory.one());
  }

  from(value: Polynomial): ExtendedFiniteField {
    return new ExtendedFiniteField(this, this.mod, value);
  }
}

export class ExtendedFiniteField implements Field<ExtendedFiniteField, Polynomial> {
  public readonly factory: ExtendedFiniteFieldFactory;
  public readonly mod: Polynomial;
  public readonly n: Polynomial;

  constructor(factory: ExtendedFiniteFieldFactory, mod: Polynomial, n: Polynomial) {
    this.factory = factory;
    this.mod = mod;
    this.n = n.mod(this.mod);
  }

  eq(other: ExtendedFiniteField): boolean {
    return this.n.eq(other.n);
  }

  clone(): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n);
  }

  isZero(): boolean {
    return this.n.isZero();
  }

  add(other: ExtendedFiniteField): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.add(other.n));
  }

  sub(other: ExtendedFiniteField): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.sub(other.n));
  }

  mul(other: ExtendedFiniteField): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.mul(other.n));
  }

  div(other: ExtendedFiniteField): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.div(other.n));
  }

  scale(n: bigint): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.scale(n));
  }

  pow(n: bigint): ExtendedFiniteField {
    return new ExtendedFiniteField(this.factory, this.mod, this.n.pow(n));
  }
}
