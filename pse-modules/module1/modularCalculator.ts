type Operation = "+" | "-" | "*";

function modularCalculator(op: Operation, num1: number, num2: number, mod: number) {
  if (op === "+") return (num1 + num2) % mod;
  if (op === "-") return (num1 - num2 + mod) % mod; // + mod to avoid negative numbers
  if (op === "*") return (num1 * num2) % mod;
}

if (import.meta.vitest) {
  it("modularCalculator", async () => {
    expect(modularCalculator("+", 10, 15, 12)).toEqual(1);
    expect(modularCalculator("-", 10, 15, 12)).toEqual(7);
    expect(modularCalculator("*", 10, 15, 12)).toEqual(6);
  });
}
