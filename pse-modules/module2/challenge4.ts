class PedersenCommitment {
  p: bigint;
  g: bigint;
  h: bigint | null;
  r: bigint | null;
  s: bigint | null;

  constructor() {
    // Set prime number (p) and generator (g)
    this.p = BigInt(23); // use a large prime in a real-world scenario
    this.g = BigInt(4); // use a large number in a real-world scenario
    this.h = null;
    this.r = null;
    this.s = null;
  }

  // Generate 'h' with a random number 'r' (h = g^r mod p)
  generateH() {
    // TODO: Generate a random number r (and save it to this.r)
    // TODO: Calculate h using g, r and p (and save it to this.h)
    this.r = BigInt(Math.floor(Math.random() * Number(this.p)));
    this.h = this.g ** this.r % this.p; // use fast modular exponentiation in a real-world scenario
  }

  // Generate the commitment (g^s * h^r mod p)
  generateCommitment(s: number) {
    // TODO: Convert s to BigInt (and save it to this.s)
    // TODO: Calculate and return the commitment using g, s, h, r and p
    this.s = BigInt(s);
    return (this.g ** this.s * (this.h as bigint) ** (this.r as bigint)) % this.p;
  }

  // Reveal the secret number and random number (s, r)
  reveal() {
    // TODO: Return the secret and random number
    return { s: this.s as bigint, r: this.r as bigint };
  }

  // Verify the commitment (g^s * h^r mod p)
  verify(s: bigint, r: bigint, C: bigint) {
    const reval = (this.g ** s * (this.h as bigint) ** r) % this.p;
    return reval === C;
  }
}

// Test the PedersenCommitment
const pc = new PedersenCommitment();
pc.generateH();

// Party A: Generate a commitment
let secretNumber = 7;
let commitment = pc.generateCommitment(secretNumber);
console.log("Commitment: ", commitment);

// Party A: Reveal the secret and random number
let reveal = pc.reveal();
console.log("Revealed: ", reveal);

// Party B: Verify the commitment
let verification = pc.verify(reveal.s, reveal.r, commitment);
console.log("Verification: ", verification);
