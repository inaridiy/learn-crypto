/**
 * 一旦暗号学的に安全ではないランダム関数を実装する
 * 0以上p未満の乱数を返す
 * @param p
 */
export const random = (p: bigint) => {
  const randomHexes = Array(5)
    .fill("0")
    .map(() => Math.random().toString(16).slice(2))
    .join("");
  return BigInt("0x" + randomHexes) % p;
};
