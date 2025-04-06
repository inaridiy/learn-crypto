import { Field, FieldFactory } from "./interface";
import { extendedGCDBigint, modBigint, pow } from "../math";

export class FiniteFieldFactory implements FieldFactory<FiniteField, bigint> {
  constructor(public readonly p: bigint) {}

  zero(): FiniteField {
    return new FiniteField(this, this.p, 0n);
  }

  one(): FiniteField {
    return new FiniteField(this, this.p, 1n);
  }

  from(value: bigint): FiniteField {
    return new FiniteField(this, this.p, value);
  }
}

export class FiniteField implements Field<FiniteField, bigint> {
  public readonly factory: FiniteFieldFactory;
  public readonly p: bigint;
  public readonly n: bigint;

  constructor(factory: FiniteFieldFactory, p: bigint, n: bigint) {
    this.factory = factory;
    this.p = p;
    this.n = modBigint(n, this.p);
  }

  eq(other: FiniteField): boolean {
    return this.n === other.n;
  }

  clone(): FiniteField {
    return new FiniteField(this.factory, this.p, this.n);
  }

  isZero(): boolean {
    return this.n === 0n;
  }

  add(other: FiniteField): FiniteField {
    return new FiniteField(this.factory, this.p, this.n + other.n);
  }

  sub(other: FiniteField): FiniteField {
    return new FiniteField(this.factory, this.p, this.n - other.n);
  }

  mul(other: FiniteField): FiniteField {
    return new FiniteField(this.factory, this.p, this.n * other.n);
  }

  scale(n: bigint): FiniteField {
    return new FiniteField(this.factory, this.p, this.n * n);
  }

  inverse(): FiniteField {
    const [gcd, x] = extendedGCDBigint(this.n, this.p);
    if (gcd !== 1n) throw new Error("No inverse");
    return new FiniteField(this.factory, this.p, x);
  }

  div(other: FiniteField): FiniteField {
    if (other.n === 0n) throw new Error("Division by zero");
    const otherInv = other.inverse();
    return this.mul(otherInv);
  }

  pow(n: bigint): FiniteField {
    return pow<FiniteField>(this, n);
  }
}
