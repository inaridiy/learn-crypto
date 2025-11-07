import { FiniteField, FiniteFieldElement } from "./finite-field";
import { PolynomialFactory, PolynomialOnFF } from "./polynomial";

export const lagrangeInterpolationFF = (
  factory: PolynomialFactory<FiniteFieldElement, bigint>,
  points: (readonly [FiniteFieldElement, FiniteFieldElement])[]
) => {
  if (points.length === 0) return factory.zero();

  const field = factory.coeffField;
  const xPoly = factory.from([field.zero(), field.one()]);

  const basisPolynomials = points.map(([xi, yi], i) => {
    const others = points.filter((_, j) => j !== i);

    const numerator = others
      .map(([xj]) => xPoly.sub(factory.from([xj])))
      .reduce((acc, poly) => acc.mul(poly), factory.one());

    const denominator = others
      .map(([xj]) => xi.sub(xj))
      .reduce((acc, diff) => acc.mul(diff), field.one());

    const scaledBasis = numerator.mul(factory.from([yi.div(denominator)]));

    return scaledBasis;
  });

  return basisPolynomials.reduce((acc, poly) => acc.add(poly), factory.zero());
};

export const matrixToPolynomials = (
  field: FiniteField,
  matrix: FiniteFieldElement[][]
): PolynomialOnFF[] => {
  if (matrix.length === 0) return [];

  const width = matrix[0].length;
  const factory = new PolynomialFactory(field);
  const evaluationPoints = matrix.map((_, idx) => field.from(BigInt(idx + 1)));

  return Array.from({ length: width }, (_, column) => {
    const samples = matrix.map((row, rowIndex) => {
      const value = row[column] ?? field.zero();
      return [evaluationPoints[rowIndex], value] as const;
    });

    return lagrangeInterpolationFF(factory, samples);
  });
};

if (import.meta.vitest) {
  const { describe, test, expect } = import.meta.vitest;
  const F = new FiniteField(13n);
  const PF = new PolynomialFactory(F);

  describe("lagrangeInterpolationFF", () => {
    test("returns zero polynomial for empty input", () => {
      const result = lagrangeInterpolationFF(PF, []);
      expect(result.isZero()).toBe(true);
    });

    test("reconstructs polynomial from sample points", () => {
      const poly = PF.from([1n, 2n, 3n]);
      const points = [0n, 1n, 3n].map((value) => {
        const x = F.from(value);
        return [x, poly.eval(x)] as [FiniteFieldElement, FiniteFieldElement];
      });

      const result = lagrangeInterpolationFF(PF, points);

      expect(result.eq(poly)).toBe(true);
    });

    test("interpolated polynomial matches evaluations", () => {
      const points: [FiniteFieldElement, FiniteFieldElement][] = [
        [F.from(2n), F.from(5n)],
        [F.from(5n), F.from(8n)],
        [F.from(7n), F.from(4n)],
      ];

      const poly = lagrangeInterpolationFF(PF, points);

      points.forEach(([x, y]) => {
        expect(poly.eval(x).eq(y)).toBe(true);
      });
    });
  });
}
