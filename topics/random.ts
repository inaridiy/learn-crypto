/**
 * 一旦暗号学的に安全ではないランダム関数を実装する
 * 0以上p未満の乱数を返す
 * @param p
 */
export const random = (p: bigint) => {
  return BigInt((Math.random() * Number(p)) | 0);
};
