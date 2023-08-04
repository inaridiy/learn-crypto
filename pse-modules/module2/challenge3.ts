import { MerkleTree } from "merkletreejs";
import crypto from "crypto";

function hashFunction(data: string) {
  return crypto.createHash("sha256").update(data).digest();
}

// Create tree
const leaves = ["a", "b", "c", "d"].map((x) => hashFunction(x));
// TODO: Build the Merkle tree using the leaves and hashFunction. Compute the root of the tree and print it.
const tree = new MerkleTree(leaves, hashFunction);
const root = tree.getRoot().toString("hex");
console.log(root);

// Generate and verify proof
const leaf = hashFunction("b");
// TODO: Generate a proof for the leaf 'b' and verify it against the root of the tree. It should return true if the leaf is part of the tree.
const proof = tree.getProof(leaf);
const verified = tree.verify(proof, leaf, root);
console.log(verified);
