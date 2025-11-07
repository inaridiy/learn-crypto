import { FiniteFieldElement } from "../primitive/finite-field";
import { PolynomialFactory, PolynomialOnFF } from "../primitive/polynomial";
import { R1CS, WitnessLike } from "./r1cs";

const buildLagrangeBasis = (
  domain: FiniteFieldElement[],
  factory: PolynomialFactory<FiniteFieldElement, bigint>
) => {
  const basis: PolynomialOnFF[] = [];

  for (let i = 0; i < domain.length; i++) {
    let numerator = factory.one();
    let denominator = domain[i].structure.one();

    for (let j = 0; j < domain.length; j++) {
      if (i === j) continue;
      const xMinusTj = factory.from([domain[j].negate(), domain[j].structure.one()]);
      numerator = numerator.mul(xMinusTj);
      denominator = denominator.mul(domain[i].sub(domain[j]));
    }

    const denomInv = denominator.inverse();
    const scalar = factory.from([denomInv]);
    basis.push(numerator.mul(scalar));
  }

  return basis;
};

const buildWirePolynomials = (
  matrix: FiniteFieldElement[][],
  basis: PolynomialOnFF[],
  factory: PolynomialFactory<FiniteFieldElement, bigint>
) => {
  const columns = matrix[0]?.length ?? 0;
  const polynomials: PolynomialOnFF[] = [];
  for (let col = 0; col < columns; col++) {
    let poly = factory.zero();
    for (let row = 0; row < matrix.length; row++) {
      const coeff = matrix[row][col];
      if (!coeff || coeff.isZero()) continue;
      const scaledBasis = basis[row].mul(factory.from([coeff]));
      poly = poly.add(scaledBasis);
    }
    polynomials.push(poly);
  }
  return polynomials;
};

const buildTargetPolynomial = (
  domain: FiniteFieldElement[],
  factory: PolynomialFactory<FiniteFieldElement, bigint>
) =>
  domain.reduce(
    (acc, point) => acc.mul(factory.from([point.negate(), point.structure.one()])),
    factory.one()
  );

export class QAP {
  public readonly factory: PolynomialFactory<FiniteFieldElement, bigint>;
  public readonly domain: FiniteFieldElement[];
  public readonly targetPolynomial: PolynomialOnFF;
  public readonly leftPolynomials: PolynomialOnFF[];
  public readonly rightPolynomials: PolynomialOnFF[];
  public readonly outputPolynomials: PolynomialOnFF[];

  private constructor(public readonly r1cs: R1CS) {
    if (r1cs.constraintCount() === 0) throw new Error("R1CS must have at least one constraint");

    this.factory = new PolynomialFactory(r1cs.field);
    this.domain = Array.from({ length: r1cs.constraintCount() }, (_, idx) =>
      r1cs.field.from(BigInt(idx + 1))
    );

    const basis = buildLagrangeBasis(this.domain, this.factory);
    this.leftPolynomials = buildWirePolynomials(r1cs.A, basis, this.factory);
    this.rightPolynomials = buildWirePolynomials(r1cs.B, basis, this.factory);
    this.outputPolynomials = buildWirePolynomials(r1cs.C, basis, this.factory);
    this.targetPolynomial = buildTargetPolynomial(this.domain, this.factory);
  }

  static fromR1CS(r1cs: R1CS) {
    return new QAP(r1cs);
  }

  witnessPolynomials(witness: WitnessLike) {
    const w = this.r1cs.coerceWitness(witness);
    return {
      A: this.combine(this.leftPolynomials, w),
      B: this.combine(this.rightPolynomials, w),
      C: this.combine(this.outputPolynomials, w),
    } as const;
  }

  verifyWitness(witness: WitnessLike) {
    const { A, B, C } = this.witnessPolynomials(witness);
    const numerator = A.mul(B).sub(C);
    return numerator.remainder(this.targetPolynomial).isZero();
  }

  evaluateAtTau(witness: WitnessLike, tau: FiniteFieldElement) {
    const { A, B, C } = this.witnessPolynomials(witness);
    return {
      a: A.eval(tau),
      b: B.eval(tau),
      c: C.eval(tau),
      t: this.targetPolynomial.eval(tau),
    } as const;
  }

  private combine(polynomials: PolynomialOnFF[], witness: FiniteFieldElement[]) {
    let acc = this.factory.zero();
    for (let i = 0; i < polynomials.length; i++) {
      const coeff = witness[i];
      if (!coeff || coeff.isZero()) continue;
      const scaled = polynomials[i].mul(this.factory.from([coeff]));
      acc = acc.add(scaled);
    }
    return acc;
  }
}
