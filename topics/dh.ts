/**
 * https://zenn.dev/herumi/articles/sd202203-ecc-1
 */

import { print } from "../utils/print";
import { random } from "./random";

//事前に決めた素数
const p = 65537n;
const g = 3n;

const a = random(p); //Aliceの秘密
const b = random(p); //Bobの秘密

print("State", { p, g, a, b });

const A = g ** a % p; //AliceがBobに渡すやつ
const B = g ** b % p; //BobがAliceに渡すやつ

print("Public", { A, B });

const s1 = B ** a % p; //AliceがBobからもらったBを使って計算する
const s2 = A ** b % p; //Bobが計算する

print("ShardSecret", { s1, s2, OK: s1 === s2 });
