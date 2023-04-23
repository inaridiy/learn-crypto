import { print } from "../../utils/print";
import { EllipticCurve } from "./EllipticCurve";
import { random } from "./../random";

const BIG_NUM = 2n ** 256n;

const curve = EllipticCurve.SECP256K1; //標準的な楕円曲線,a=0,b=7
const G = EllipticCurve.SECP256K1_G; //基点
const a = random(BIG_NUM); //Aliceの秘密鍵
const b = random(BIG_NUM); //Bobの秘密鍵

print("Private", { a, b });

const aP = curve.multiply(G, a); //Aliceの公開鍵
const bP = curve.multiply(G, b); //Bobの公開鍵

print("Public", { aP, bP });

const baP = curve.multiply(aP, b); //共有鍵 aP * b = bP * a (交換法則)
const abP = curve.multiply(bP, a); //共有鍵

print("Shared", { baP, abP, OK: baP.x === abP.x });
