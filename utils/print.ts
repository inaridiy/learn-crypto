import { Polynomial } from "../topics/paring/curve/Polynomial";

const SPLITTER_LENGTH = 64;

const printableTransform = (x: any) => {
  if (typeof x === "bigint") return "0x" + x.toString(16);
  if (typeof x === "object" && typeof x.serializeToHexStr === "function")
    return "0x" + x.serializeToHexStr().slice(0, 20) + "...";
  if (x instanceof Polynomial) return x.toString();
  return x;
};

const printRecordNested = (title: string, x: Record<string, any>, depth: number) => {
  if (title.length === 0) console.log("=".repeat(SPLITTER_LENGTH));
  else console.log("== " + title + " " + "=".repeat(SPLITTER_LENGTH - title.length - 4));

  for (const key in x) {
    if (x[key] instanceof Polynomial) console.log(`${key}: ${printableTransform(x[key])}`);
    else if (typeof x[key] === "object" && !x[key].serializeToHexStr) {
      for (const key2 in x[key]) {
        console.log(`${key}.${key2}: ${printableTransform(x[key][key2])}`);
      }
    } else {
      console.log(`${key}: ${printableTransform(x[key])}`);
    }
  }
  console.log("=".repeat(SPLITTER_LENGTH));
};

type print = (title: any, x?: Record<string, any>) => void;
export const print: print = (title, x) => {
  if (typeof title === "string" && x) printRecordNested(title, x, 0);
  else printRecordNested("", title, 0);
};
