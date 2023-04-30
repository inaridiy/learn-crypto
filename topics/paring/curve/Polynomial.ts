export class Polynomial {
  public readonly coefficients: number[];

  static normalize(coefficients: number[], degree = coefficients.length): number[] {
    let result = [...coefficients];
    while (result.length > 1 && result[result.length - 1] === 0) result.pop();
    return result;
  }

  static mustBePolynomial(polynomial: Polynomial | number[]): Polynomial {
    if (polynomial instanceof Polynomial) return polynomial;
    else if (Array.isArray(polynomial)) return new Polynomial(polynomial);
    else throw new Error("Not a polynomial");
  }

  constructor(coefficients_: number[]) {
    this.coefficients = Polynomial.normalize(coefficients_);
  }

  add(other_: Polynomial | number[]): Polynomial {
    const other = Polynomial.mustBePolynomial(other_);
    const maxLen = Math.max(this.coefficients.length, other.coefficients.length);
    const result = Array(maxLen).fill(0);

    for (let i = 0; i < maxLen; i++) {
      const a = this.coefficients[i] || 0;
      const b = other.coefficients[i] || 0;
      result[i] = a + b;
    }

    return new Polynomial(result);
  }

  sub(other_: Polynomial | number[]): Polynomial {
    const other = Polynomial.mustBePolynomial(other_);
    const maxLen = Math.max(this.coefficients.length, other.coefficients.length);
    const result = Array(maxLen).fill(0);

    for (let i = 0; i < maxLen; i++) {
      const a = this.coefficients[i] || 0;
      const b = other.coefficients[i] || 0;
      result[i] = a - b;
    }

    return new Polynomial(result);
  }

  mul(other_: Polynomial | number[]): Polynomial {
    const other = Polynomial.mustBePolynomial(other_);
    const resultDegree = this.coefficients.length + other.coefficients.length - 1;
    const result = Array(resultDegree).fill(0);

    for (let i = 0; i < this.coefficients.length; i++) {
      for (let j = 0; j < other.coefficients.length; j++) {
        result[i + j] += this.coefficients[i] * other.coefficients[j];
      }
    }

    return new Polynomial(result);
  }

  eval(x: number): number {
    let result = 0;
    for (let i = this.coefficients.length - 1; i >= 0; i--) {
      result = result * x + this.coefficients[i];
    }
    return result;
  }

  toString(): string {
    const terms = this.coefficients.map((c, i) => {
      if (c === 0) return "";
      else if (i === 0) return `${c}`;
      else if (i === 1) return `${c}x`;
      else return `${c}x^${i}`;
    });
    return terms.reverse().join(" + ");
  }
}
