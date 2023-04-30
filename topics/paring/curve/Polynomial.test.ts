import { Polynomial } from "./Polynomial";

describe("Polynomial", () => {
  it("add Polynomial", () => {
    const p1 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
    const p2 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2

    expect(p1.add(p2).coefficients).toEqual([2, 4, 6]); // 2 + 4x + 6x^2
    expect(p1.add([4, 5, 6]).coefficients).toEqual([5, 7, 9]); // 5 + 7x + 9x^2

    const p3 = new Polynomial([1, 2, 3, 4]); // 1 + 2x + 3x^2 + 4x^3
    const p4 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2

    expect(p3.add(p4).coefficients).toEqual([2, 4, 6, 4]); // 2 + 4x + 6x^2 + 4x^3

    const p5 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
    const p6 = new Polynomial([-1, 2, -3, -4]); // -1 + 2x - 3x^2 - 4x^3

    expect(p5.add(p6).coefficients).toEqual([0, 4, 0, -4]); // 0 + 4x + 0x^2 + 0x^3
  });

  it("sub Polynomial", () => {
    const p1 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
    const p2 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2

    expect(p1.sub(p2).coefficients).toEqual([0]); // 0 + 0x + 0x^2
    expect(p1.sub([4, 5, 6]).coefficients).toEqual([-3, -3, -3]); // -3 + -3x + -3x^2

    const p3 = new Polynomial([1, 2, 3, 4]); // 1 + 2x + 3x^2 + 4x^3
    const p4 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2

    expect(p3.sub(p4).coefficients).toEqual([0, 0, 0, 4]); // 0 + 0x + 0x^2 + 4x^3
  });

  it("mul Polynomial", () => {
    const p1 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
    const p2 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2

    expect(p1.mul(p2).coefficients).toEqual([1, 4, 10, 12, 9]); // 1 + 4x + 10x^2 + 12x^3 + 9x^4
    expect(p1.mul([4, 5, 6]).coefficients).toEqual([4, 13, 28, 27, 18]); // 4 + 13x + 28x^2 + 27x^3 + 18x^4

    const p3 = new Polynomial([1, 2, 3, 4]); // 1 + 2x + 3x^2 + 4x^3
    const p4 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2

    expect(p3.mul(p4).coefficients).toEqual([1, 4, 10, 16, 17, 12]); // 1 + 4x + 10x^2 + 16x^3 + 17x^4 + 12x^5

    const p5 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
    const p6 = new Polynomial([1, 2, 3, 4]); // 1 + 2x + 3x^2 + 4x^3

    expect(p5.mul(p6).coefficients).toEqual([1, 4, 10, 16, 17, 12]); // 1 + 4x + 10x^2 + 16x^3 + 17x^4 + 12x^5

    const p7 = new Polynomial([-1, -2, -3]);
    const p8 = new Polynomial([1, 2, 3]);

    expect(p7.mul(p8).coefficients).toEqual([-1, -4, -10, -12, -9]);
  });

  it("eval Polynomial", () => {
    const p1 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
    const x1 = 1;

    expect(p1.eval(x1)).toEqual(6); // 1 + 2 + 3

    const p2 = new Polynomial([1, 2, 3, 4]); // 1 + 2x + 3x^2 + 4x^3
    const x2 = 2;

    expect(p2.eval(x2)).toEqual(49); // 1 + 4 + 12 + 32
  });

  it("toString", () => {
    const p1 = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
    expect(p1.toString()).toEqual("3x^2 + 2x + 1");

    const p2 = new Polynomial([1, 2, 3, 4]); // 1 + 2x + 3x^2 + 4x^3
    expect(p2.toString()).toEqual("4x^3 + 3x^2 + 2x + 1");
  });
});
