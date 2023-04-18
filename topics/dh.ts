/**
 * https://zenn.dev/herumi/articles/sd202203-ecc-1
 * 
このコードは、Diffie-Hellman鍵共有 (DH暗号) の簡易的な実装を示しています。このアルゴリズムでは、AliceとBobは互いに秘密鍵と公開鍵を交換し、最終的な共有シークレットを計算します。s1とs2が一致する理由は、数学的な性質に基づいています。

まず、アルゴリズムにおいて、以下のような計算が行われています。

Aliceの秘密鍵（a）と公開鍵（A）: A = g^a mod p
Bobの秘密鍵（b）と公開鍵（B）: B = g^b mod p
AliceとBobは、互いに公開鍵を交換します。

Aliceは、Bobの公開鍵（B）を使って共有シークレット（s1）を計算します: s1 = B^a mod p
Bobは、Aliceの公開鍵（A）を使って共有シークレット（s2）を計算します: s2 = A^b mod p
ここで、s1とs2が一致する理由を示します。

s1 = B^a mod p = (g^b mod p)^a mod p = g^(ab) mod p
s2 = A^b mod p = (g^a mod p)^b mod p = g^(ab) mod p

s1とs2はどちらも g^(ab) mod p を計算しているため、一致します。この性質により、AliceとBobは互いに同じ共有シークレットを計算できます。これがDiffie-Hellman鍵共有アルゴリズムの基本原理です。
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
