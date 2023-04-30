import { FQ } from "./FQ";
import { Field } from "../../simple-ec/EllipticCurve";

describe("FQ", () => {
  it("add FQ", () => {
    const p = 11n;
    const a = new FQ(p, 3n);
    const b = new FQ(p, 5n);
    const c = a.add(b);
    expect(c.n).toBe(8n);

    const d = a.add(5n);
    expect(d.n).toBe(8n);

    const e = a.add(9n);
    expect(e.n).toBe(1n);
  });

  it("mul FQ", () => {
    const p = 11n;
    const a = new FQ(p, 3n);
    const b = new FQ(p, 5n);
    const c = a.mul(b);
    expect(c.n).toBe(4n);

    const d = a.mul(5n);
    expect(d.n).toBe(4n);

    const e = a.mul(9n);
    expect(e.n).toBe(5n);
  });

  it("div FQ", () => {
    const sample = new Field(11n);

    const p = 11n;
    const a = new FQ(p, 3n);
    const b = new FQ(p, 5n);
    const c = a.div(b);

    expect(c.n).toBe(sample.div(3n, 5n));
  });
});
