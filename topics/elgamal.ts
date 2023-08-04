import { print } from "../utils/print";
import { random } from "./random";
import { fastPow, modInverse } from "./math";

const p = 2n ** 127n - 1n; //素数
const g = 2n; //生成元

const encrypt = (sP: bigint, m: bigint) => {
  const t = random(p - 4n) + 2n;
  const c1 = fastPow(g, t, p); //暗号文(1
  const c2 = (fastPow(sP, t, p) * m) % p; //暗号文(2
  return [c1, c2] as const;
};

const decrypt = (s: bigint, c1: bigint, c2: bigint) => {
  const dec = c2 * modInverse(fastPow(c1, s, p), p);
  return (dec + p) % p;
};

const s = random(p - 4n) + 2n; //秘密鍵
const sP = fastPow(g, s, p); //公開鍵
const msg = 0x123456789n;

print("State", { s, sP, msg });

const [c1, c2] = encrypt(sP, msg);
print("Encrypted", { c1, c2 });

const dec = decrypt(s, c1, c2);
print("Decrypted", { dec, Ok: dec === msg });
