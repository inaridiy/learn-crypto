import { bench, run, group } from "mitata";
import { FQ } from "./FQ";

const BIG_PRIME = 2n ** 256n - 2n ** 32n - 977n;
const A = new FQ(BIG_PRIME, 3n).div(5n);
const SMALL = 3n;
const MEDIAN = 2n ** 3n;
const BIG = 2n ** 8n;

group("FQ binary pow", () => {
  bench("small", () => {
    A.pow(SMALL);
  });
  bench("median", () => {
    A.pow(MEDIAN);
  });
  bench("big", () => {
    A.pow(BIG);
  });
});

group("FQ just pow", () => {
  bench("small", () => {
    A.mod(A.n ** SMALL);
  });
  bench("median", () => {
    A.mod(A.n ** MEDIAN);
  });
  bench("big", () => {
    A.mod(A.n ** BIG);
  });
});

run();
