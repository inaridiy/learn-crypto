//! Pedersen commitment に対する logarithmic dot-product argument。
//!
//! Hyrax 論文 Appendix A.3 (Figures 7--8) の `prooflog-of-dot-prod` を、
//! Fiat--Shamir 変換で非対話化したもの。公開ベクトル `a` と、commit された
//! ベクトル `x` / スカラー `y` に対して `y = <x, a>` を zero knowledge で示す。

use ark_ec::CurveGroup;
use ark_ff::{Field, Zero};
use ark_std::{
    UniformRand,
    rand::{CryptoRng, Rng},
};

use crate::primitive::{ScalarPedersen, Transcript, VectorPedersen, fold_halves, inner_product};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InnerProductProof<G: CurveGroup> {
    pub minus: Vec<G>,
    pub plus: Vec<G>,

    pub delta: G,
    pub beta: G,

    pub z1: G::ScalarField,
    pub z2: G::ScalarField,
}

fn append_statement<G: CurveGroup>(
    transcript: &mut Transcript,
    vector_committer: &VectorPedersen<G>,
    scalar_committer: &ScalarPedersen<G>,
    public_vector: &[G::ScalarField],
    vector_commitment: &G,
    scalar_commitment: &G,
) {
    transcript.append_usize(b"ipa-vector-length", public_vector.len());
    for generator in &vector_committer.generators {
        transcript.append_serializable(b"ipa-vector-generator", generator);
    }
    transcript.append_serializable(b"ipa-scalar-generator", &scalar_committer.generator);
    transcript.append_serializable(b"ipa-blind-generator", &vector_committer.blind);
    for value in public_vector {
        transcript.append_serializable(b"ipa-public-vector", value);
    }
    transcript.append_serializable(b"ipa-vector-commitment", vector_commitment);
    transcript.append_serializable(b"ipa-scalar-commitment", scalar_commitment);
}

fn challenge_nonzero<G: CurveGroup>(transcript: &mut Transcript) -> G::ScalarField {
    // Figure 7 samples c from the non-zero field elements. A zero Fiat--Shamir
    // output is committed to the transcript and deterministically retried.
    loop {
        let challenge = transcript.challenge_field::<G::ScalarField>(b"ipa-challenge");
        if !challenge.is_zero() {
            return challenge;
        }
    }
}

impl<G: CurveGroup> InnerProductProof<G> {
    /// `y = <x, a>` に対する Figure 7--8 の logarithmic proof を作る。
    #[allow(clippy::too_many_arguments)]
    pub fn prove(
        vector_committer: &VectorPedersen<G>,
        scalar_committer: &ScalarPedersen<G>,
        vector_commitment: &G,
        scalar_commitment: &G,
        x: &[G::ScalarField],
        vector_blind: &G::ScalarField,
        a: &[G::ScalarField],
        scalar_blind: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> Self {
        Self::assert_shape(vector_committer, scalar_committer, x, a);

        let y = inner_product(x, a);
        assert_eq!(
            vector_committer.commit(x, vector_blind),
            *vector_commitment,
            "vector commitment does not match its opening"
        );
        assert_eq!(
            scalar_committer.commit(&y, scalar_blind),
            *scalar_commitment,
            "scalar commitment does not commit to the inner product"
        );

        append_statement(
            transcript,
            vector_committer,
            scalar_committer,
            a,
            vector_commitment,
            scalar_commitment,
        );

        let mut x = x.to_vec();
        let mut a = a.to_vec();
        let mut generators = vector_committer.generators.clone();
        let mut combined_blind = *vector_blind + scalar_blind;
        let mut minus = Vec::with_capacity(x.len().ilog2() as usize);
        let mut plus = Vec::with_capacity(x.len().ilog2() as usize);

        while x.len() > 1 {
            let half = x.len() / 2;
            let (x_l, x_r) = x.split_at(half);
            let (a_l, a_r) = a.split_at(half);
            let (g_l, g_r) = generators.split_at(half);

            let blind_minus = G::ScalarField::rand(rng);
            let blind_plus = G::ScalarField::rand(rng);
            let m_minus = inner_product(g_r, x_l)
                + scalar_committer.generator * inner_product(x_l, a_r)
                + vector_committer.blind * blind_minus;
            let m_plus = inner_product(g_l, x_r)
                + scalar_committer.generator * inner_product(x_r, a_l)
                + vector_committer.blind * blind_plus;

            transcript.append_serializable(b"ipa-minus", &m_minus);
            transcript.append_serializable(b"ipa-plus", &m_plus);
            minus.push(m_minus);
            plus.push(m_plus);

            let c = challenge_nonzero::<G>(transcript);
            let c_inv = c.inverse().expect("non-zero challenge has an inverse");

            combined_blind =
                c.square() * blind_minus + combined_blind + c_inv.square() * blind_plus;
            x = fold_halves(&x, c, c_inv);
            a = fold_halves(&a, c_inv, c);
            generators = fold_halves(&generators, c_inv, c);
        }

        let d = G::ScalarField::rand(rng);
        let blind_delta = G::ScalarField::rand(rng);
        let blind_beta = G::ScalarField::rand(rng);
        let delta = generators[0] * d + vector_committer.blind * blind_delta;
        let beta = scalar_committer.generator * d + vector_committer.blind * blind_beta;
        transcript.append_serializable(b"ipa-delta", &delta);
        transcript.append_serializable(b"ipa-beta", &beta);

        let c = challenge_nonzero::<G>(transcript);
        let y_final = x[0] * a[0];
        let z1 = d + c * y_final;
        let z2 = a[0] * (c * combined_blind + blind_beta) + blind_delta;

        Self {
            minus,
            plus,
            delta,
            beta,
            z1,
            z2,
        }
    }

    /// Figure 7 の reduction と Figure 8 の Equation (14) を検証する。
    pub fn verify(
        &self,
        vector_committer: &VectorPedersen<G>,
        scalar_committer: &ScalarPedersen<G>,
        public_vector: &[G::ScalarField],
        vector_commitment: &G,
        scalar_commitment: &G,
        transcript: &mut Transcript,
    ) -> bool {
        if !Self::valid_shape(vector_committer, scalar_committer, public_vector)
            || self.minus.len() != public_vector.len().ilog2() as usize
            || self.plus.len() != self.minus.len()
        {
            return false;
        }

        append_statement(
            transcript,
            vector_committer,
            scalar_committer,
            public_vector,
            vector_commitment,
            scalar_commitment,
        );

        let mut statement = *vector_commitment + scalar_commitment;
        let mut a = public_vector.to_vec();
        let mut generators = vector_committer.generators.clone();

        for (m_minus, m_plus) in self.minus.iter().zip(&self.plus) {
            transcript.append_serializable(b"ipa-minus", m_minus);
            transcript.append_serializable(b"ipa-plus", m_plus);

            let c = challenge_nonzero::<G>(transcript);
            let c_inv = c.inverse().expect("non-zero challenge has an inverse");

            statement = *m_minus * c.square() + statement + *m_plus * c_inv.square();
            a = fold_halves(&a, c_inv, c);
            generators = fold_halves(&generators, c_inv, c);
        }

        transcript.append_serializable(b"ipa-delta", &self.delta);
        transcript.append_serializable(b"ipa-beta", &self.beta);
        let c = challenge_nonzero::<G>(transcript);

        // Additive notation for Equation (14):
        // a_hat (c Upsilon_hat + beta) + delta
        //   = z1 (g_hat + a_hat g) + z2 h.
        (statement * c + self.beta) * a[0] + self.delta
            == (generators[0] + scalar_committer.generator * a[0]) * self.z1
                + vector_committer.blind * self.z2
    }

    fn assert_shape(
        vector_committer: &VectorPedersen<G>,
        scalar_committer: &ScalarPedersen<G>,
        x: &[G::ScalarField],
        a: &[G::ScalarField],
    ) {
        assert!(
            Self::valid_shape(vector_committer, scalar_committer, a) && x.len() == a.len(),
            "IPA vectors must have the same non-zero power-of-two length as the generators, and share h"
        );
    }

    fn valid_shape(
        vector_committer: &VectorPedersen<G>,
        scalar_committer: &ScalarPedersen<G>,
        a: &[G::ScalarField],
    ) -> bool {
        !a.is_empty()
            && a.len().is_power_of_two()
            && vector_committer.len() == a.len()
            && vector_committer.blind == scalar_committer.blind
    }
}

#[cfg(test)]
mod tests {
    use super::InnerProductProof;
    use crate::primitive::{ScalarPedersen, Transcript, VectorPedersen, inner_product};
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use ark_std::rand::{SeedableRng, rngs::StdRng};
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    fn test_rng() -> StdRng {
        StdRng::seed_from_u64(42)
    }

    fn committers(n: usize) -> (VectorPedersen<G1Projective>, ScalarPedersen<G1Projective>) {
        (
            VectorPedersen::setup::<G1Hasher>(b"ipa-test", n).unwrap(),
            ScalarPedersen::setup::<G1Hasher>(b"ipa-test").unwrap(),
        )
    }

    fn prove_example(
        n: usize,
    ) -> (
        VectorPedersen<G1Projective>,
        ScalarPedersen<G1Projective>,
        Vec<F>,
        G1Projective,
        G1Projective,
        InnerProductProof<G1Projective>,
    ) {
        let (vector, scalar) = committers(n);
        let x = (1..=n).map(|i| F::from(i as u64)).collect::<Vec<_>>();
        let a = (0..n)
            .map(|i| F::from((2 * i + 3) as u64))
            .collect::<Vec<_>>();
        let vector_blind = F::from(11);
        let scalar_blind = F::from(13);
        let y = inner_product(&x, &a);
        let c_x = vector.commit(&x, &vector_blind);
        let c_y = scalar.commit(&y, &scalar_blind);
        let mut transcript = Transcript::new(b"ipa-proof");
        let proof = InnerProductProof::prove(
            &vector,
            &scalar,
            &c_x,
            &c_y,
            &x,
            &vector_blind,
            &a,
            &scalar_blind,
            &mut transcript,
            &mut test_rng(),
        );
        (vector, scalar, a, c_x, c_y, proof)
    }

    #[test]
    fn logarithmic_ipa_accepts_valid_statements_including_length_one() {
        for n in [1, 2, 4, 8] {
            let (vector, scalar, a, c_x, c_y, proof) = prove_example(n);
            let mut transcript = Transcript::new(b"ipa-proof");
            assert!(proof.verify(&vector, &scalar, &a, &c_x, &c_y, &mut transcript));
            assert_eq!(proof.minus.len(), n.ilog2() as usize);
        }
    }

    #[test]
    fn logarithmic_ipa_rejects_changed_public_inputs_and_proof() {
        let (vector, scalar, mut a, c_x, c_y, proof) = prove_example(4);
        a[0] += F::from(1);
        let mut transcript = Transcript::new(b"ipa-proof");
        assert!(!proof.verify(&vector, &scalar, &a, &c_x, &c_y, &mut transcript));

        let (vector, scalar, a, c_x, c_y, mut proof) = prove_example(4);
        let wrong_c_y = c_y + scalar.generator;
        let mut transcript = Transcript::new(b"ipa-proof");
        assert!(!proof.verify(&vector, &scalar, &a, &c_x, &wrong_c_y, &mut transcript));

        proof.z1 += F::from(1);
        let mut transcript = Transcript::new(b"ipa-proof");
        assert!(!proof.verify(&vector, &scalar, &a, &c_x, &c_y, &mut transcript));
    }

    #[test]
    fn verifier_rejects_wrong_shapes_and_distinct_blind_generators() {
        let (vector, scalar, a, c_x, c_y, mut proof) = prove_example(4);
        proof.minus.pop();
        let mut transcript = Transcript::new(b"ipa-proof");
        assert!(!proof.verify(&vector, &scalar, &a, &c_x, &c_y, &mut transcript));

        let other_scalar = ScalarPedersen::setup::<G1Hasher>(b"other-domain").unwrap();
        let (_, _, _, _, _, proof) = prove_example(4);
        let mut transcript = Transcript::new(b"ipa-proof");
        assert!(!proof.verify(&vector, &other_scalar, &a, &c_x, &c_y, &mut transcript));
    }
}
