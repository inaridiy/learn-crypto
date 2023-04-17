const SPLITTER_LENGTH = 64;

const printRecordNested = (title: string, x: Record<string, any>, depth: number) => {
  if (title.length === 0) console.log("=".repeat(SPLITTER_LENGTH));
  else console.log("== " + title + " " + "=".repeat(SPLITTER_LENGTH - title.length - 4));

  for (const key in x) {
    if (typeof x[key] === "object") {
      printRecordNested(key, x[key], depth + 1);
    } else {
      console.log(`${key}: ${x[key]}`);
    }
  }
  console.log("=".repeat(SPLITTER_LENGTH));
};

type print = (title: any, x?: Record<string, any>) => void;
export const print: print = (title, x) => {
  if (typeof title === "string" && x) printRecordNested(title, x, 0);
  else printRecordNested("", title, 0);
};
