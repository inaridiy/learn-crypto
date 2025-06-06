# inaridiy's Cryptography Learning Notes

## Overview

Repository where inaridiy implements various protocols related to cryptography. No library!
※I also implement PSE assignments in this repository and use the library there as well.

## Execution

Testing Cryptographic Protocols

```bash
npm run test
```

Some things need to be run directly.

```
npx tsx path/to/file.ts
```

## Implemented Protocols

### 1. Diffie-Hellman Key Exchange

[Simple DH Key Exchange](./topics/dh.ts)  
[ECDH Key Exchange](./topics/simple-ec/ecdh.ts)

### 2. Elgamal Encryption

[Simple Elgamal Encryption](./topics/elgamal.ts)  
[EC Elgamal Encryption](./topics/simple-ec/ecelgamal.ts)  
[Lifted Elgamal Encryption and additive quasi-homomorphism](./topics/simple-ec/liftedecelgamal.ts)

### 3. DSA/ECDSA Digital Signature

[ECDSA Digital Signature](./topics/simple-ec/ecdsa.ts)

### 4. Paring-Based Cryptography

[Check bilinear](./topics/paring/bilinear.ts)  
[BLS Signature](./topics/paring/bls.ts)  
[KZG Commitment](./topics/kzg/kzg.ts)

### 5. Elliptic Curve Cryptography

[Elliptic Curve Arithmetic](./topics/simple-ec/EllipticCurve/EllipticCurve.ts)  
[ECDH Key Exchange](./topics/simple-ec/ecdh.ts)  
[EC Elgamal Encryption](./topics/simple-ec/ecelgamal.ts)  
[Lifted Elgamal Encryption and additive quasi-homomorphism](./topics/simple-ec/liftedecelgamal.ts)  
[ECDSA Digital Signature](./topics/simple-ec/ecdsa.ts)

## PSE Summer Contribution Program Exercise

[Module 1](./pse-modules/module1)
[Module 2](./pse-modules/module2/)
