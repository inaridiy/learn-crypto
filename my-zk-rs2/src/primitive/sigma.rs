//! Pedersen commitment 上の Σ-protocol 群 (Fiat-Shamir で非対話化済み)。
//!
//! すべて $\mathrm{Com}(v; r) = v G + r H$ (スカラー) /
//! $\mathrm{Com}(\vec{v}; r) = \sum_i v_i G_i + r H$ (ベクトル) という
//! Pedersen commitment に対する honest-verifier zero-knowledge な証明で、
//! zk sum-check や SpartanNIZK の最終検証をコミットメントの上だけで行うために使う。

use ark_ec::CurveGroup;
use ark_std::{UniformRand, rand::Rng};

use super::helpers::inner_product;
use super::pedersen::Pedersen;
use super::transcript::Transcript;

/// Commitment $C = x G + r H$ の開示値 $(x, r)$ を知っていることの証明
/// (zero-knowledge proof of knowledge, Schnorr protocol)。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KnowledgeProof<G: CurveGroup> {
    pub alpha: G,
    pub z1: G::ScalarField,
    pub z2: G::ScalarField,
}

impl<G: CurveGroup> KnowledgeProof<G> {
    pub fn prove(
        scalar_committer: &Pedersen<G>,
        x: &G::ScalarField,
        blind_x: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut impl Rng,
    ) -> Self {
        let c_x = scalar_committer.commit_scalar(x, blind_x);
        transcript.append_serializable(b"knowledge-c-x", &c_x);

        let (t1, t2) = (G::ScalarField::rand(rng), G::ScalarField::rand(rng));
        let alpha = scalar_committer.commit_scalar(&t1, &t2);
        transcript.append_serializable(b"knowledge-alpha", &alpha);

        let c: G::ScalarField = transcript.challenge_field(b"knowledge-challenge");

        Self {
            alpha,
            z1: c * x + t1,
            z2: c * blind_x + t2,
        }
    }

    /// $z_1 G + z_2 H = c \cdot C + \alpha$ を確認する。
    pub fn verify(
        &self,
        scalar_committer: &Pedersen<G>,
        c_x: &G,
        transcript: &mut Transcript,
    ) -> bool {
        transcript.append_serializable(b"knowledge-c-x", c_x);
        transcript.append_serializable(b"knowledge-alpha", &self.alpha);

        let c: G::ScalarField = transcript.challenge_field(b"knowledge-challenge");

        scalar_committer.commit_scalar(&self.z1, &self.z2) == *c_x * c + self.alpha
    }
}

/// $C_1$ と $C_2$ が同じ値への commitment であることの証明。
///
/// $C_1 - C_2 = (r_1 - r_2) H$ なので、$C_1 - C_2$ の $H$ を底とした
/// 離散対数 (blind の差) の知識を Schnorr protocol で証明すればよい。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EqualityProof<G: CurveGroup> {
    pub alpha: G,
    pub z: G::ScalarField,
}

impl<G: CurveGroup> EqualityProof<G> {
    /// `blind_diff` $= r_1 - r_2$ を知る prover が証明を作る。
    pub fn prove(
        blind_generator: &G,
        c1: &G,
        c2: &G,
        blind_diff: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut impl Rng,
    ) -> Self {
        transcript.append_serializable(b"equality-c1", c1);
        transcript.append_serializable(b"equality-c2", c2);

        let r = G::ScalarField::rand(rng);
        let alpha = *blind_generator * r;
        transcript.append_serializable(b"equality-alpha", &alpha);

        let c: G::ScalarField = transcript.challenge_field(b"equality-challenge");
        let z = c * blind_diff + r;

        Self { alpha, z }
    }

    /// $z H = \alpha + c (C_1 - C_2)$ を確認する。
    pub fn verify(&self, blind_generator: &G, c1: &G, c2: &G, transcript: &mut Transcript) -> bool {
        transcript.append_serializable(b"equality-c1", c1);
        transcript.append_serializable(b"equality-c2", c2);
        transcript.append_serializable(b"equality-alpha", &self.alpha);

        let c: G::ScalarField = transcript.challenge_field(b"equality-challenge");

        *blind_generator * self.z == self.alpha + (*c1 - c2) * c
    }
}

/// $C_x, C_y, C_z$ に対して $z = x \cdot y$ であることの証明 (product argument)。
///
/// $C_x = x G + r_x H,\; C_y = y G + r_y H,\; C_z = x y G + r_z H$。
/// 鍵となる観察は $C_z = y \cdot C_x + (r_z - y r_x) H$ と書けること。
/// つまり「$C_y$ の値 $y$」と「$C_x$ を底とした表現」の両方に同じ $y$ が
/// 現れることを 3 本の Schnorr 型の等式で同時に証明する。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProductProof<G: CurveGroup> {
    pub alpha: G,
    pub beta: G,
    pub delta: G,
    pub z: [G::ScalarField; 5],
}

impl<G: CurveGroup> ProductProof<G> {
    #[allow(clippy::too_many_arguments)]
    pub fn prove(
        scalar_committer: &Pedersen<G>,
        x: &G::ScalarField,
        blind_x: &G::ScalarField,
        y: &G::ScalarField,
        blind_y: &G::ScalarField,
        blind_z: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut impl Rng,
    ) -> Self {
        let c_x = scalar_committer.commit_scalar(x, blind_x);
        let c_y = scalar_committer.commit_scalar(y, blind_y);
        let c_z = scalar_committer.commit_scalar(&(*x * y), blind_z);
        transcript.append_serializable(b"product-c-x", &c_x);
        transcript.append_serializable(b"product-c-y", &c_y);
        transcript.append_serializable(b"product-c-z", &c_z);

        let (b1, b2, b3, b4, b5) = (
            G::ScalarField::rand(rng),
            G::ScalarField::rand(rng),
            G::ScalarField::rand(rng),
            G::ScalarField::rand(rng),
            G::ScalarField::rand(rng),
        );

        // alpha, beta は C_x, C_y への Schnorr commitment、
        // delta は「C_x を底にした C_z の表現」への Schnorr commitment。
        let alpha = scalar_committer.commit_scalar(&b1, &b2);
        let beta = scalar_committer.commit_scalar(&b3, &b4);
        let delta = c_x * b3 + scalar_committer.blind * b5;
        transcript.append_serializable(b"product-alpha", &alpha);
        transcript.append_serializable(b"product-beta", &beta);
        transcript.append_serializable(b"product-delta", &delta);

        let c: G::ScalarField = transcript.challenge_field(b"product-challenge");

        let z = [
            b1 + c * x,
            b2 + c * blind_x,
            b3 + c * y,
            b4 + c * blind_y,
            b5 + c * (*blind_z - *blind_x * y),
        ];

        Self {
            alpha,
            beta,
            delta,
            z,
        }
    }

    pub fn verify(
        &self,
        scalar_committer: &Pedersen<G>,
        c_x: &G,
        c_y: &G,
        c_z: &G,
        transcript: &mut Transcript,
    ) -> bool {
        transcript.append_serializable(b"product-c-x", c_x);
        transcript.append_serializable(b"product-c-y", c_y);
        transcript.append_serializable(b"product-c-z", c_z);
        transcript.append_serializable(b"product-alpha", &self.alpha);
        transcript.append_serializable(b"product-beta", &self.beta);
        transcript.append_serializable(b"product-delta", &self.delta);

        let c: G::ScalarField = transcript.challenge_field(b"product-challenge");
        let [z1, z2, z3, z4, z5] = self.z;

        scalar_committer.commit_scalar(&z1, &z2) == self.alpha + *c_x * c
            && scalar_committer.commit_scalar(&z3, &z4) == self.beta + *c_y * c
            && *c_x * z3 + scalar_committer.blind * z5 == self.delta + *c_z * c
    }
}

/// commit 済みベクトル $\vec{x}$ と公開ベクトル $\vec{a}$ の内積が
/// commit 済みスカラー $y$ に等しいこと、つまり
/// $y = \langle \vec{x}, \vec{a} \rangle$ の証明。
///
/// - $C_x = \mathrm{Com}(\vec{x}; r_x)$ (vector committer)
/// - $C_y = \mathrm{Com}(y; r_y)$ (scalar committer)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DotProductProof<G: CurveGroup> {
    pub delta: G,
    pub beta: G,
    pub z: Vec<G::ScalarField>,
    pub z_delta: G::ScalarField,
    pub z_beta: G::ScalarField,
}

impl<G: CurveGroup> DotProductProof<G> {
    #[allow(clippy::too_many_arguments)]
    pub fn prove(
        vec_committer: &Pedersen<G>,
        scalar_committer: &Pedersen<G>,
        x: &[G::ScalarField],
        blind_x: &G::ScalarField,
        a: &[G::ScalarField],
        blind_y: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut impl Rng,
    ) -> Self {
        let y = inner_product(x, a);
        let c_x = vec_committer.commit(x, blind_x);
        let c_y = scalar_committer.commit_scalar(&y, blind_y);
        Self::append_statement(transcript, a, &c_x, &c_y);

        // ランダムベクトル d への commitment (delta) と <d, a> への commitment (beta)。
        let d: Vec<G::ScalarField> = (0..x.len()).map(|_| G::ScalarField::rand(rng)).collect();
        let blind_delta = G::ScalarField::rand(rng);
        let blind_beta = G::ScalarField::rand(rng);
        let delta = vec_committer.commit(&d, &blind_delta);
        let beta = scalar_committer.commit_scalar(&inner_product(&d, a), &blind_beta);
        transcript.append_serializable(b"dot-product-delta", &delta);
        transcript.append_serializable(b"dot-product-beta", &beta);

        let c: G::ScalarField = transcript.challenge_field(b"dot-product-challenge");

        let z = x.iter().zip(&d).map(|(x_i, d_i)| c * x_i + d_i).collect();
        let z_delta = c * blind_x + blind_delta;
        let z_beta = c * blind_y + blind_beta;

        Self {
            delta,
            beta,
            z,
            z_delta,
            z_beta,
        }
    }

    /// $\mathrm{Com}(\vec{z}; z_\delta) = c \cdot C_x + \delta$ と
    /// $\mathrm{Com}(\langle \vec{z}, \vec{a} \rangle; z_\beta) = c \cdot C_y + \beta$
    /// を確認する。
    pub fn verify(
        &self,
        vec_committer: &Pedersen<G>,
        scalar_committer: &Pedersen<G>,
        a: &[G::ScalarField],
        c_x: &G,
        c_y: &G,
        transcript: &mut Transcript,
    ) -> bool {
        if self.z.len() != a.len() {
            return false;
        }

        Self::append_statement(transcript, a, c_x, c_y);
        transcript.append_serializable(b"dot-product-delta", &self.delta);
        transcript.append_serializable(b"dot-product-beta", &self.beta);

        let c: G::ScalarField = transcript.challenge_field(b"dot-product-challenge");

        vec_committer.commit(&self.z, &self.z_delta) == *c_x * c + self.delta
            && scalar_committer.commit_scalar(&inner_product(&self.z, a), &self.z_beta)
                == *c_y * c + self.beta
    }

    fn append_statement(transcript: &mut Transcript, a: &[G::ScalarField], c_x: &G, c_y: &G) {
        for value in a {
            transcript.append_field(b"dot-product-a", value);
        }
        transcript.append_serializable(b"dot-product-c-x", c_x);
        transcript.append_serializable(b"dot-product-c-y", c_y);
    }
}

#[cfg(test)]
mod tests {
    use super::{DotProductProof, EqualityProof, KnowledgeProof, ProductProof};
    use crate::primitive::{Pedersen, Transcript, inner_product};
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_std::test_rng;
    use sha2::Sha256;

    type G1Hasher = MapToCurveBasedHasher<
        G1Projective,
        ark_ff::field_hashers::DefaultFieldHasher<Sha256, 128>,
        WBMap<Config>,
    >;

    fn committers() -> (Pedersen<G1Projective>, Pedersen<G1Projective>) {
        (
            Pedersen::setup::<G1Hasher>(b"sigma-test", "vec", 4).unwrap(),
            Pedersen::setup::<G1Hasher>(b"sigma-test", "scalar", 1).unwrap(),
        )
    }

    #[test]
    fn knowledge_proof_accepts_a_known_opening() {
        let (_, scalar) = committers();
        let mut rng = test_rng();
        let (value, blind) = (F::from(42), F::from(7));
        let c = scalar.commit_scalar(&value, &blind);

        let mut prover_transcript = Transcript::new(b"sigma-test");
        let proof =
            KnowledgeProof::prove(&scalar, &value, &blind, &mut prover_transcript, &mut rng);

        let mut verifier_transcript = Transcript::new(b"sigma-test");
        assert!(proof.verify(&scalar, &c, &mut verifier_transcript));
    }

    #[test]
    fn knowledge_proof_rejects_a_different_commitment() {
        let (_, scalar) = committers();
        let mut rng = test_rng();
        let (value, blind) = (F::from(42), F::from(7));
        let other_c = scalar.commit_scalar(&F::from(43), &blind);

        let mut prover_transcript = Transcript::new(b"sigma-test");
        let proof =
            KnowledgeProof::prove(&scalar, &value, &blind, &mut prover_transcript, &mut rng);

        let mut verifier_transcript = Transcript::new(b"sigma-test");
        assert!(!proof.verify(&scalar, &other_c, &mut verifier_transcript));
    }

    #[test]
    fn equality_proof_accepts_commitments_to_the_same_value() {
        let (_, scalar) = committers();
        let mut rng = test_rng();
        let (value, r1, r2) = (F::from(42), F::from(3), F::from(5));
        let c1 = scalar.commit_scalar(&value, &r1);
        let c2 = scalar.commit_scalar(&value, &r2);

        let mut prover_transcript = Transcript::new(b"sigma-test");
        let proof = EqualityProof::prove(
            &scalar.blind,
            &c1,
            &c2,
            &(r1 - r2),
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"sigma-test");
        assert!(proof.verify(&scalar.blind, &c1, &c2, &mut verifier_transcript));
    }

    #[test]
    fn equality_proof_rejects_commitments_to_different_values() {
        let (_, scalar) = committers();
        let mut rng = test_rng();
        let (r1, r2) = (F::from(3), F::from(5));
        let c1 = scalar.commit_scalar(&F::from(42), &r1);
        let c2 = scalar.commit_scalar(&F::from(43), &r2);

        let mut prover_transcript = Transcript::new(b"sigma-test");
        let proof = EqualityProof::prove(
            &scalar.blind,
            &c1,
            &c2,
            &(r1 - r2),
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"sigma-test");
        assert!(!proof.verify(&scalar.blind, &c1, &c2, &mut verifier_transcript));
    }

    #[test]
    fn product_proof_accepts_a_true_product() {
        let (_, scalar) = committers();
        let mut rng = test_rng();
        let (x, y) = (F::from(6), F::from(7));
        let (r_x, r_y, r_z) = (F::from(11), F::from(13), F::from(17));
        let c_x = scalar.commit_scalar(&x, &r_x);
        let c_y = scalar.commit_scalar(&y, &r_y);
        let c_z = scalar.commit_scalar(&(x * y), &r_z);

        let mut prover_transcript = Transcript::new(b"sigma-test");
        let proof = ProductProof::prove(
            &scalar,
            &x,
            &r_x,
            &y,
            &r_y,
            &r_z,
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"sigma-test");
        assert!(proof.verify(&scalar, &c_x, &c_y, &c_z, &mut verifier_transcript));
    }

    #[test]
    fn product_proof_rejects_a_wrong_product() {
        let (_, scalar) = committers();
        let mut rng = test_rng();
        let (x, y) = (F::from(6), F::from(7));
        let (r_x, r_y, r_z) = (F::from(11), F::from(13), F::from(17));
        let c_x = scalar.commit_scalar(&x, &r_x);
        let c_y = scalar.commit_scalar(&y, &r_y);
        let wrong_c_z = scalar.commit_scalar(&(x * y + F::from(1)), &r_z);

        let mut prover_transcript = Transcript::new(b"sigma-test");
        let proof = ProductProof::prove(
            &scalar,
            &x,
            &r_x,
            &y,
            &r_y,
            &r_z,
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"sigma-test");
        assert!(!proof.verify(&scalar, &c_x, &c_y, &wrong_c_z, &mut verifier_transcript));
    }

    #[test]
    fn dot_product_proof_accepts_a_true_inner_product() {
        let (vec, scalar) = committers();
        let mut rng = test_rng();
        let x = [F::from(1), F::from(2), F::from(3), F::from(4)];
        let a = [F::from(5), F::from(6), F::from(7), F::from(8)];
        let (r_x, r_y) = (F::from(11), F::from(13));
        let c_x = vec.commit(&x, &r_x);
        let c_y = scalar.commit_scalar(&inner_product(&x, &a), &r_y);

        let mut prover_transcript = Transcript::new(b"sigma-test");
        let proof = DotProductProof::prove(
            &vec,
            &scalar,
            &x,
            &r_x,
            &a,
            &r_y,
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"sigma-test");
        assert!(proof.verify(&vec, &scalar, &a, &c_x, &c_y, &mut verifier_transcript));
    }

    #[test]
    fn dot_product_proof_rejects_a_wrong_inner_product() {
        let (vec, scalar) = committers();
        let mut rng = test_rng();
        let x = [F::from(1), F::from(2), F::from(3), F::from(4)];
        let a = [F::from(5), F::from(6), F::from(7), F::from(8)];
        let (r_x, r_y) = (F::from(11), F::from(13));
        let c_x = vec.commit(&x, &r_x);
        let wrong_c_y = scalar.commit_scalar(&(inner_product(&x, &a) + F::from(1)), &r_y);

        let mut prover_transcript = Transcript::new(b"sigma-test");
        let proof = DotProductProof::prove(
            &vec,
            &scalar,
            &x,
            &r_x,
            &a,
            &r_y,
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"sigma-test");
        assert!(!proof.verify(
            &vec,
            &scalar,
            &a,
            &c_x,
            &wrong_c_y,
            &mut verifier_transcript
        ));
    }
}
