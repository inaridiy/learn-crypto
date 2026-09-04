//! Pedersen commitment に対する sigma protocol 群。すべて Fiat--Shamir 変換で非対話化する。
//!
//! - [`InnerProductProof`][]: Hyrax 論文 Appendix A.3 (Figures 7--8) の
//!   `prooflog-of-dot-prod`。公開ベクトル `a` と、commit されたベクトル `x` /
//!   スカラー `y` に対して `y = <x, a>` を対数サイズで示す。
//! - [`LinearInnerProductProof`][]: 同じ関係を線形サイズで示す単純な版。
//!   ZK sumcheck のラウンド多項式のように `x` が短いときに使う。
//! - [`KnowledgeProof`][]: commitment の opening を知っていることの Schnorr 証明。
//! - [`EqualityProof`][]: 二つの commitment が同じ値を隠していることの証明。
//! - [`ProductProof`][]: $z = x y$ を三つの commitment 上で示す証明。
//!
//! 後半の 4 つは Spartan の参照実装(`nizk` モジュール)に対応する。

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
    transcript.append_serializable(b"ipa-vector-blind-generator", &vector_committer.blind);
    transcript.append_serializable(b"ipa-scalar-blind-generator", &scalar_committer.blind);
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

/// 線形サイズの dot-product argument。
///
/// `y = <x, a>` を、`x` の commitment $C_x$ と `y` の commitment $C_y$ に対して示す。
/// prover は乱数ベクトル $\vec{d}$ を commit した $\delta$ と $\langle \vec{a}, \vec{d} \rangle$ を
/// commit した $\beta$ を送り、challenge $c$ に対して $\vec{z} = c \vec{x} + \vec{d}$ を開く。
/// 証明サイズは $O(|x|)$ なので、`x` が短い(sumcheck のラウンド多項式の係数など)場合に使う。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearInnerProductProof<G: CurveGroup> {
    pub delta: G,
    pub beta: G,
    pub z: Vec<G::ScalarField>,
    pub z_delta: G::ScalarField,
    pub z_beta: G::ScalarField,
}

impl<G: CurveGroup> LinearInnerProductProof<G> {
    /// `y = <x, a>` に対する線形サイズの proof を作る。引数は [`InnerProductProof::prove`] と同じ。
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
        assert!(
            Self::valid_shape(vector_committer, a) && x.len() == a.len(),
            "linear IPA vectors must have the same non-zero length as the generators"
        );
        assert_eq!(
            vector_committer.commit(x, vector_blind),
            *vector_commitment,
            "vector commitment does not match its opening"
        );
        assert_eq!(
            scalar_committer.commit(&inner_product(x, a), scalar_blind),
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

        let d = (0..x.len())
            .map(|_| G::ScalarField::rand(rng))
            .collect::<Vec<_>>();
        let r_delta = G::ScalarField::rand(rng);
        let r_beta = G::ScalarField::rand(rng);
        let delta = vector_committer.commit(&d, &r_delta);
        let beta = scalar_committer.commit(&inner_product(&d, a), &r_beta);
        transcript.append_serializable(b"lipa-delta", &delta);
        transcript.append_serializable(b"lipa-beta", &beta);

        let c = transcript.challenge_field::<G::ScalarField>(b"lipa-challenge");
        let z = x.iter().zip(&d).map(|(x, d)| c * x + d).collect();

        Self {
            delta,
            beta,
            z,
            z_delta: c * vector_blind + r_delta,
            z_beta: c * scalar_blind + r_beta,
        }
    }

    /// $c C_x + \delta = \mathrm{Com}(\vec{z}; z_\delta)$ と
    /// $c C_y + \beta = \mathrm{Com}(\langle \vec{z}, \vec{a} \rangle; z_\beta)$ を検証する。
    pub fn verify(
        &self,
        vector_committer: &VectorPedersen<G>,
        scalar_committer: &ScalarPedersen<G>,
        public_vector: &[G::ScalarField],
        vector_commitment: &G,
        scalar_commitment: &G,
        transcript: &mut Transcript,
    ) -> bool {
        if !Self::valid_shape(vector_committer, public_vector)
            || self.z.len() != public_vector.len()
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
        transcript.append_serializable(b"lipa-delta", &self.delta);
        transcript.append_serializable(b"lipa-beta", &self.beta);
        let c = transcript.challenge_field::<G::ScalarField>(b"lipa-challenge");

        vector_committer.commit(&self.z, &self.z_delta) == *vector_commitment * c + self.delta
            && scalar_committer.commit(&inner_product(&self.z, public_vector), &self.z_beta)
                == *scalar_commitment * c + self.beta
    }

    fn valid_shape(vector_committer: &VectorPedersen<G>, a: &[G::ScalarField]) -> bool {
        !a.is_empty() && vector_committer.len() == a.len()
    }
}

/// Scalar commitment の鍵と、それに対する statement の commitment 列を transcript に bind する。
fn append_scalar_statement<G: CurveGroup>(
    transcript: &mut Transcript,
    protocol: &[u8],
    committer: &ScalarPedersen<G>,
    commitments: &[&G],
) {
    transcript.append_bytes(b"sigma-protocol", protocol);
    transcript.append_serializable(b"sigma-generator", &committer.generator);
    transcript.append_serializable(b"sigma-blind-generator", &committer.blind);
    for commitment in commitments {
        transcript.append_serializable(b"sigma-commitment", *commitment);
    }
}

/// Pedersen commitment $C = v G + r H$ の opening $(v, r)$ を知っていることの Schnorr 証明。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KnowledgeProof<G: CurveGroup> {
    pub alpha: G,
    pub z1: G::ScalarField,
    pub z2: G::ScalarField,
}

impl<G: CurveGroup> KnowledgeProof<G> {
    pub fn prove(
        committer: &ScalarPedersen<G>,
        commitment: &G,
        value: &G::ScalarField,
        blind: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> Self {
        assert_eq!(
            committer.commit(value, blind),
            *commitment,
            "commitment does not match its opening"
        );

        let t1 = G::ScalarField::rand(rng);
        let t2 = G::ScalarField::rand(rng);
        let alpha = committer.commit(&t1, &t2);

        append_scalar_statement(transcript, b"knowledge", committer, &[commitment]);
        transcript.append_serializable(b"knowledge-alpha", &alpha);
        let c = transcript.challenge_field::<G::ScalarField>(b"knowledge-challenge");

        Self {
            alpha,
            z1: t1 + c * value,
            z2: t2 + c * blind,
        }
    }

    /// $\mathrm{Com}(z_1; z_2) = c C + \alpha$ を検証する。
    pub fn verify(
        &self,
        committer: &ScalarPedersen<G>,
        commitment: &G,
        transcript: &mut Transcript,
    ) -> bool {
        append_scalar_statement(transcript, b"knowledge", committer, &[commitment]);
        transcript.append_serializable(b"knowledge-alpha", &self.alpha);
        let c = transcript.challenge_field::<G::ScalarField>(b"knowledge-challenge");

        committer.commit(&self.z1, &self.z2) == *commitment * c + self.alpha
    }
}

/// 二つの Pedersen commitment $C_1, C_2$ が同じ値を隠していることの証明。
///
/// $C_1 - C_2 = (r_1 - r_2) H$ なので、$H$ についての離散対数 $r_1 - r_2$ の知識証明になる。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EqualityProof<G: CurveGroup> {
    pub alpha: G,
    pub z: G::ScalarField,
}

impl<G: CurveGroup> EqualityProof<G> {
    #[allow(clippy::too_many_arguments)]
    pub fn prove(
        committer: &ScalarPedersen<G>,
        lhs_commitment: &G,
        rhs_commitment: &G,
        value: &G::ScalarField,
        lhs_blind: &G::ScalarField,
        rhs_blind: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> Self {
        assert_eq!(
            committer.commit(value, lhs_blind),
            *lhs_commitment,
            "left commitment does not match its opening"
        );
        assert_eq!(
            committer.commit(value, rhs_blind),
            *rhs_commitment,
            "right commitment does not match its opening"
        );

        let t = G::ScalarField::rand(rng);
        let alpha = committer.blind * t;

        append_scalar_statement(
            transcript,
            b"equality",
            committer,
            &[lhs_commitment, rhs_commitment],
        );
        transcript.append_serializable(b"equality-alpha", &alpha);
        let c = transcript.challenge_field::<G::ScalarField>(b"equality-challenge");

        Self {
            alpha,
            z: t + c * (*lhs_blind - rhs_blind),
        }
    }

    /// $z H = c (C_1 - C_2) + \alpha$ を検証する。
    pub fn verify(
        &self,
        committer: &ScalarPedersen<G>,
        lhs_commitment: &G,
        rhs_commitment: &G,
        transcript: &mut Transcript,
    ) -> bool {
        append_scalar_statement(
            transcript,
            b"equality",
            committer,
            &[lhs_commitment, rhs_commitment],
        );
        transcript.append_serializable(b"equality-alpha", &self.alpha);
        let c = transcript.challenge_field::<G::ScalarField>(b"equality-challenge");

        committer.blind * self.z == (*lhs_commitment - rhs_commitment) * c + self.alpha
    }
}

/// $X = \mathrm{Com}(x)$, $Y = \mathrm{Com}(y)$, $Z = \mathrm{Com}(z)$ に対して $z = x y$ を示す証明。
///
/// $X$ と $Y$ の opening の知識証明に加え、$Z$ を「$X$ を generator とみなした
/// $y$ の commitment」$Z = y X + (r_Z - r_X y) H$ として開くことで積の関係を結び付ける。
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
        committer: &ScalarPedersen<G>,
        x_commitment: &G,
        y_commitment: &G,
        z_commitment: &G,
        x: &G::ScalarField,
        x_blind: &G::ScalarField,
        y: &G::ScalarField,
        y_blind: &G::ScalarField,
        z_blind: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> Self {
        assert_eq!(
            committer.commit(x, x_blind),
            *x_commitment,
            "x commitment does not match its opening"
        );
        assert_eq!(
            committer.commit(y, y_blind),
            *y_commitment,
            "y commitment does not match its opening"
        );
        assert_eq!(
            committer.commit(&(*x * y), z_blind),
            *z_commitment,
            "z commitment does not commit to the product"
        );

        let b: [G::ScalarField; 5] = core::array::from_fn(|_| G::ScalarField::rand(rng));
        let alpha = committer.commit(&b[0], &b[1]);
        let beta = committer.commit(&b[2], &b[3]);
        let delta = *x_commitment * b[2] + committer.blind * b[4];

        append_scalar_statement(
            transcript,
            b"product",
            committer,
            &[x_commitment, y_commitment, z_commitment],
        );
        transcript.append_serializable(b"product-alpha", &alpha);
        transcript.append_serializable(b"product-beta", &beta);
        transcript.append_serializable(b"product-delta", &delta);
        let c = transcript.challenge_field::<G::ScalarField>(b"product-challenge");

        Self {
            alpha,
            beta,
            delta,
            z: [
                b[0] + c * x,
                b[1] + c * x_blind,
                b[2] + c * y,
                b[3] + c * y_blind,
                b[4] + c * (*z_blind - *x_blind * y),
            ],
        }
    }

    /// $\mathrm{Com}(z_1; z_2) = c X + \alpha$、$\mathrm{Com}(z_3; z_4) = c Y + \beta$、
    /// $z_3 X + z_5 H = c Z + \delta$ を検証する。
    pub fn verify(
        &self,
        committer: &ScalarPedersen<G>,
        x_commitment: &G,
        y_commitment: &G,
        z_commitment: &G,
        transcript: &mut Transcript,
    ) -> bool {
        append_scalar_statement(
            transcript,
            b"product",
            committer,
            &[x_commitment, y_commitment, z_commitment],
        );
        transcript.append_serializable(b"product-alpha", &self.alpha);
        transcript.append_serializable(b"product-beta", &self.beta);
        transcript.append_serializable(b"product-delta", &self.delta);
        let c = transcript.challenge_field::<G::ScalarField>(b"product-challenge");

        let [z1, z2, z3, z4, z5] = self.z;
        committer.commit(&z1, &z2) == *x_commitment * c + self.alpha
            && committer.commit(&z3, &z4) == *y_commitment * c + self.beta
            && *x_commitment * z3 + committer.blind * z5 == *z_commitment * c + self.delta
    }
}

#[cfg(test)]
mod tests {
    use super::{
        EqualityProof, InnerProductProof, KnowledgeProof, LinearInnerProductProof, ProductProof,
    };
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
    #[test]
    fn linear_ipa_accepts_valid_statements_and_rejects_tampering() {
        for n in [1, 3, 4] {
            let (vector, scalar) = committers(n);
            let x = (1..=n).map(|i| F::from(i as u64)).collect::<Vec<_>>();
            let a = (0..n).map(|i| F::from((i + 2) as u64)).collect::<Vec<_>>();
            let (vector_blind, scalar_blind) = (F::from(11), F::from(13));
            let c_x = vector.commit(&x, &vector_blind);
            let c_y = scalar.commit(&inner_product(&x, &a), &scalar_blind);
            let mut transcript = Transcript::new(b"lipa-proof");
            let proof = LinearInnerProductProof::prove(
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
            assert_eq!(proof.z.len(), n);

            let mut transcript = Transcript::new(b"lipa-proof");
            assert!(proof.verify(&vector, &scalar, &a, &c_x, &c_y, &mut transcript));

            let mut wrong_a = a.clone();
            wrong_a[0] += F::from(1);
            let mut transcript = Transcript::new(b"lipa-proof");
            assert!(!proof.verify(&vector, &scalar, &wrong_a, &c_x, &c_y, &mut transcript));

            let mut transcript = Transcript::new(b"lipa-proof");
            assert!(!proof.verify(
                &vector,
                &scalar,
                &a,
                &(c_x + vector.generators[0]),
                &c_y,
                &mut transcript
            ));

            let mut tampered = proof.clone();
            tampered.z_beta += F::from(1);
            let mut transcript = Transcript::new(b"lipa-proof");
            assert!(!tampered.verify(&vector, &scalar, &a, &c_x, &c_y, &mut transcript));
        }
    }

    #[test]
    fn knowledge_proof_round_trips_and_binds_the_commitment() {
        let (_, scalar) = committers(1);
        let (value, blind) = (F::from(21), F::from(5));
        let commitment = scalar.commit(&value, &blind);

        let mut transcript = Transcript::new(b"knowledge");
        let proof = KnowledgeProof::prove(
            &scalar,
            &commitment,
            &value,
            &blind,
            &mut transcript,
            &mut test_rng(),
        );

        let mut transcript = Transcript::new(b"knowledge");
        assert!(proof.verify(&scalar, &commitment, &mut transcript));

        let mut transcript = Transcript::new(b"knowledge");
        assert!(!proof.verify(&scalar, &(commitment + scalar.generator), &mut transcript));

        let mut tampered = proof.clone();
        tampered.z1 += F::from(1);
        let mut transcript = Transcript::new(b"knowledge");
        assert!(!tampered.verify(&scalar, &commitment, &mut transcript));
    }

    #[test]
    fn equality_proof_accepts_same_value_and_rejects_different_values() {
        let (_, scalar) = committers(1);
        let value = F::from(21);
        let (lhs_blind, rhs_blind) = (F::from(5), F::from(7));
        let lhs = scalar.commit(&value, &lhs_blind);
        let rhs = scalar.commit(&value, &rhs_blind);

        let mut transcript = Transcript::new(b"equality");
        let proof = EqualityProof::prove(
            &scalar,
            &lhs,
            &rhs,
            &value,
            &lhs_blind,
            &rhs_blind,
            &mut transcript,
            &mut test_rng(),
        );

        let mut transcript = Transcript::new(b"equality");
        assert!(proof.verify(&scalar, &lhs, &rhs, &mut transcript));

        // 同じ blind 差でも値が異なると、C_1 - C_2 が H の倍数でなくなる。
        let other = scalar.commit(&(value + F::from(1)), &rhs_blind);
        let mut transcript = Transcript::new(b"equality");
        assert!(!proof.verify(&scalar, &lhs, &other, &mut transcript));
    }

    #[test]
    fn product_proof_accepts_products_and_rejects_non_products() {
        let (_, scalar) = committers(1);
        let (x, y) = (F::from(6), F::from(7));
        let (x_blind, y_blind, z_blind) = (F::from(3), F::from(5), F::from(11));
        let x_com = scalar.commit(&x, &x_blind);
        let y_com = scalar.commit(&y, &y_blind);
        let z_com = scalar.commit(&(x * y), &z_blind);

        let mut transcript = Transcript::new(b"product");
        let proof = ProductProof::prove(
            &scalar,
            &x_com,
            &y_com,
            &z_com,
            &x,
            &x_blind,
            &y,
            &y_blind,
            &z_blind,
            &mut transcript,
            &mut test_rng(),
        );

        let mut transcript = Transcript::new(b"product");
        assert!(proof.verify(&scalar, &x_com, &y_com, &z_com, &mut transcript));

        let wrong_z = scalar.commit(&(x * y + F::from(1)), &z_blind);
        let mut transcript = Transcript::new(b"product");
        assert!(!proof.verify(&scalar, &x_com, &y_com, &wrong_z, &mut transcript));

        let mut tampered = proof.clone();
        tampered.z[4] += F::from(1);
        let mut transcript = Transcript::new(b"product");
        assert!(!tampered.verify(&scalar, &x_com, &y_com, &z_com, &mut transcript));
    }
}
