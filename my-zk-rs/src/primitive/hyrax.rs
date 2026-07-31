use crate::primitive::{BoolHyperCube, Matrix, Transcript, helpers::msm_with_bases, inner_product};

use ark_ec::{
    AffineRepr, CurveGroup,
    hashing::{HashToCurve, HashToCurveError},
};
use ark_ff::{Field, Zero};
use ark_std::{UniformRand, rand::Rng};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HyraxPCSProof<G: CurveGroup> {
    pub result_com: G,
    pub m_minus_vec: Vec<G>,
    pub m_plus_vec: Vec<G>,
    pub schnorr_com: (G, G),
    pub schnorr_res: (G::ScalarField, G::ScalarField),
}

// Multiliner多項式の変数数 = HALF_VARS_BITS * 2
// Multiliner多項式の評価点数 = 2 ^ {HALF_VARS_BITS * 2}
// 評価点数ののlog2 = 2 ^ HALF_VARS_BITS

#[derive(Clone, Debug)]
pub struct HyraxPCS<G: CurveGroup, const HALF_VARS_BITS: usize>
where
    [(); 1 << HALF_VARS_BITS]:,
{
    vec: [G; 1 << HALF_VARS_BITS], // generators for vector commitment
    vec_mul_bases: Vec<G::Affine>,
    scalar: G, // generator for scalar commitment
    blind: G,  // pedersen blinds for vector commitment and scalar commitment
}

pub type HyraxPCSCommitment<G: CurveGroup, const HALF_VARS_BITS: usize> = [G; 1 << HALF_VARS_BITS];

impl<G: CurveGroup, const HALF_VARS_BITS: usize> HyraxPCS<G, HALF_VARS_BITS>
where
    [(); 1 << HALF_VARS_BITS]:,
    [(); 1 << (HALF_VARS_BITS * 2)]:,
{
    #[inline]
    fn evaluation_index(row: usize, col: usize) -> usize {
        // Matches BoolPoint::to_index(): row bits are the low bits, col bits are the high bits.
        row + (col << HALF_VARS_BITS)
    }

    #[inline]
    fn pack_matrix(
        evaluations: &[G::ScalarField; 1 << (HALF_VARS_BITS * 2)],
    ) -> Matrix<G::ScalarField, HALF_VARS_BITS, HALF_VARS_BITS>
    where
        [(); 1 << HALF_VARS_BITS]:,
    {
        Matrix::new(std::array::from_fn(|row| {
            std::array::from_fn(|col| evaluations[Self::evaluation_index(row, col)])
        }))
    }

    fn setup_vec_generators<H: HashToCurve<G>>(
        hasher: &H,
    ) -> Result<[G; 1 << HALF_VARS_BITS], HashToCurveError> {
        let generators: [G; 1 << HALF_VARS_BITS] = (0..(1 << HALF_VARS_BITS))
            .map(|i| {
                Ok(hasher
                    .hash(format!("hyrax:vec#{i}").as_bytes())?
                    .into_group())
            })
            .collect::<Result<Vec<_>, HashToCurveError>>()?
            .try_into()
            .unwrap();

        Ok(generators)
    }

    pub fn setup<H: HashToCurve<G>>(domain: &[u8]) -> Result<Self, HashToCurveError> {
        let hasher = H::new(domain)?;

        let vec = Self::setup_vec_generators(&hasher)?;
        let vec_mul_bases = G::batch_convert_to_mul_base(&vec);

        let scalar = hasher.hash(b"hyrax:h")?.into_group();
        let blind = hasher.hash(b"hyrax:u")?.into_group();

        Ok(Self {
            vec,
            vec_mul_bases,
            scalar,
            blind,
        })
    }

    fn commit_vec(&self, values: &[G::ScalarField], r: &G::ScalarField) -> G {
        msm_with_bases::<G>(values, &self.vec_mul_bases) + self.blind * r
    }

    fn commit_scalar(&self, value: &G::ScalarField, r: &G::ScalarField) -> G {
        self.scalar * value + self.blind * r
    }

    /// 評価値 `value` に対する Pedersen commitment を計算する。
    ///
    /// Verifier が「opening proof の `result_com` は主張している評価値 `value` に
    /// 対応している」ことを確認するために使う（`result_com` そのものは検証していない
    /// 生の値なので、外部で期待するコミットメントと比較する必要がある）。
    #[inline]
    pub fn commit_value(&self, value: &G::ScalarField, blind: &G::ScalarField) -> G {
        self.commit_scalar(value, blind)
    }

    pub fn commit(
        &self,
        evaluations: &[G::ScalarField; 1 << (HALF_VARS_BITS * 2)],
        com_blinds: [G::ScalarField; 1 << HALF_VARS_BITS],
    ) -> HyraxPCSCommitment<G, HALF_VARS_BITS> {
        let matrix = Self::pack_matrix(evaluations);

        std::array::from_fn(|row| self.commit_vec(matrix.row(row), &com_blinds[row]))
    }

    fn append_statement(
        transcript: &mut Transcript,
        com: &HyraxPCSCommitment<G, HALF_VARS_BITS>,
        point: &[G::ScalarField; HALF_VARS_BITS * 2],
        ipa_statement: &G,
        result_com: &G,
    ) {
        transcript.append_usize(b"hyrax-half-vars-bits", HALF_VARS_BITS);
        for commitment in com {
            transcript.append_serializable(b"hyrax-commitment", commitment);
        }
        for coordinate in point {
            transcript.append_field(b"hyrax-point", coordinate);
        }
        transcript.append_field(b"hyrax-ipa-statement", ipa_statement);
        transcript.append_field(b"hyrax-result", result_com);
    }

    fn calc_l_r(
        point: &[G::ScalarField; HALF_VARS_BITS * 2],
    ) -> (
        [G::ScalarField; 1 << HALF_VARS_BITS],
        [G::ScalarField; 1 << HALF_VARS_BITS],
    ) {
        let row_point: [G::ScalarField; HALF_VARS_BITS] = std::array::from_fn(|i| point[i]);
        let col_point: [G::ScalarField; HALF_VARS_BITS] =
            std::array::from_fn(|i| point[HALF_VARS_BITS + i]);

        (
            BoolHyperCube::<HALF_VARS_BITS>::eq_evaluations(&row_point),
            BoolHyperCube::<HALF_VARS_BITS>::eq_evaluations(&col_point),
        )
    }

    pub fn prove_with_transcript(
        &self,
        com: &HyraxPCSCommitment<G, HALF_VARS_BITS>,
        evaluations: &[G::ScalarField; 1 << (HALF_VARS_BITS * 2)],
        com_blinds: &[G::ScalarField; 1 << HALF_VARS_BITS],
        point: &[G::ScalarField; HALF_VARS_BITS * 2],
        result_blind: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut impl Rng,
    ) -> HyraxPCSProof<G> {
        let matrix = Self::pack_matrix(evaluations);
        let (l, r) = Self::calc_l_r(point);

        // L^{\top} T
        let mut u = matrix.left_mul_vector(&l).to_vec();
        let mut q = r.to_vec();

        let result = inner_product(&u, &r);
        let mut blind = inner_product(com_blinds, &l) + *result_blind;

        let result_com = self.commit_scalar(&result, result_blind);
        let ipa_statement = msm_with_bases::<G>(&u, &self.vec_mul_bases)
            + self.scalar * result
            + self.blind * blind;

        Self::append_statement(transcript, com, point, &ipa_statement, &result_com);

        let mut g = self.vec.to_vec();
        let mut g_mul_bases = self.vec_mul_bases.to_vec();

        let mut m_minus_vec = Vec::with_capacity(HALF_VARS_BITS);
        let mut m_plus_vec = Vec::with_capacity(HALF_VARS_BITS);

        for _ in 0..HALF_VARS_BITS {
            let (r_minus, r_plus) = (G::ScalarField::rand(rng), G::ScalarField::rand(rng));

            let half = u.len() / 2;
            let (u_l, u_r) = u.split_at(half);
            let (q_l, q_r) = q.split_at(half);
            let (g_l, g_r) = g.split_at(half);
            let (g_l_mul_bases, g_r_mul_bases) = g_mul_bases.split_at(half);

            let m_minus = msm_with_bases::<G>(u_l, g_r_mul_bases)
                + self.scalar * inner_product(u_l, q_r)
                + self.blind * r_minus;
            let m_plus = msm_with_bases::<G>(u_r, g_l_mul_bases)
                + self.scalar * inner_product(u_r, q_l)
                + self.blind * r_plus;

            transcript.append_serializable(b"m_minus", &m_minus);
            transcript.append_serializable(b"m_plus", &m_plus);
            m_minus_vec.push(m_minus);
            m_plus_vec.push(m_plus);

            let c: G::ScalarField = transcript.challenge_field(b"challenge");
            let c_inv = c.inverse().expect("challenge is non-zero");

            let next_u: Vec<_> = (0..half).map(|i| c * u_l[i] + c_inv * u_r[i]).collect();
            let next_q: Vec<_> = (0..half).map(|i| c_inv * q_l[i] + c * q_r[i]).collect();

            let next_blind = c.square() * r_minus + blind + c_inv.square() * r_plus;

            let next_g = (0..half)
                .map(|i| g_l[i] * c_inv + g_r[i] * c)
                .collect::<Vec<_>>();

            u = next_u;
            q = next_q;
            blind = next_blind;
            g_mul_bases = G::batch_convert_to_mul_base(&next_g);
            g = next_g;
        }

        let (r1, r2, r3) = (
            G::ScalarField::rand(rng),
            G::ScalarField::rand(rng),
            G::ScalarField::rand(rng),
        );

        let com1 = g[0] * r1 + self.blind * r2;
        let com2 = self.scalar * r1 + self.blind * r3;

        transcript.append_serializable(b"com1", &com1);
        transcript.append_serializable(b"com2", &com2);

        let e: G::ScalarField = transcript.challenge_field(b"schnorr-challange");

        assert!(q[0] != G::ScalarField::zero(), "OMG");

        let z1 = r1 + e * u[0] * q[0];
        let z2 = q[0] * (e * blind + r3) + r2;

        HyraxPCSProof {
            result_com,
            m_minus_vec,
            m_plus_vec,
            schnorr_com: (com1, com2),
            schnorr_res: (z1, z2),
        }
    }

    pub fn verify_with_transcript(
        &self,
        com: &HyraxPCSCommitment<G, HALF_VARS_BITS>,
        point: &[G::ScalarField; HALF_VARS_BITS * 2],
        proof: &HyraxPCSProof<G>,
        transcript: &mut Transcript,
    ) -> bool {
        if proof.m_minus_vec.len() != HALF_VARS_BITS || proof.m_plus_vec.len() != HALF_VARS_BITS {
            return false;
        }

        let (l, r) = Self::calc_l_r(point);

        let c_projected = com
            .iter()
            .zip(l)
            .fold(G::zero(), |acc, (commitment, scalar)| {
                acc + *commitment * scalar
            });
        let mut statement = c_projected + proof.result_com;

        Self::append_statement(transcript, com, point, &statement, &proof.result_com);

        let mut q = r.to_vec();
        let mut g = self.vec.to_vec();

        for (m_minus, m_plus) in proof.m_minus_vec.iter().zip(&proof.m_plus_vec) {
            transcript.append_serializable(b"m_minus", m_minus);
            transcript.append_serializable(b"m_plus", m_plus);

            let c: G::ScalarField = transcript.challenge_field(b"challenge");
            let c_inv = c.inverse().expect("challenge is non-zero");

            statement = *m_minus * c.square() + statement + *m_plus * c_inv.square();

            let half = q.len() / 2;
            let (q_l, q_r) = q.split_at(half);
            let (g_l, g_r) = g.split_at(half);

            let next_q: Vec<_> = (0..half).map(|i| c_inv * q_l[i] + c * q_r[i]).collect();
            let next_g = (0..half)
                .map(|i| g_l[i] * c_inv + g_r[i] * c)
                .collect::<Vec<_>>();

            q = next_q;
            g = next_g;
        }

        let (com1, com2) = proof.schnorr_com;
        transcript.append_serializable(b"com1", &com1);
        transcript.append_serializable(b"com2", &com2);

        let e: G::ScalarField = transcript.challenge_field(b"schnorr-challange");
        let q_final = q[0];
        if q_final.is_zero() {
            return false;
        }
        let (z1, z2) = proof.schnorr_res;

        g[0] * z1 + self.scalar * (q_final * z1) + self.blind * z2
            == com1 + com2 * q_final + statement * (e * q_final)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::primitive::mle_from_evaluations;
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use ark_std::test_rng;
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    fn example_statement() -> (
        HyraxPCS<G1Projective, 2>,
        [F; 16],
        [F; 4],
        HyraxPCSCommitment<G1Projective, 2>,
        [F; 4],
        F,
    ) {
        let pcs = HyraxPCS::<G1Projective, 2>::setup::<G1Hasher>(b"hyrax-test").unwrap();
        let evaluations = std::array::from_fn(|i| F::from((i as u64 + 1).pow(2)));
        let blinds = [F::from(3), F::from(5), F::from(7), F::from(11)];
        let commitment = pcs.commit(&evaluations, blinds);
        let point = [F::from(2), F::from(3), F::from(5), F::from(7)];
        let result_blind = F::from(13);

        (pcs, evaluations, blinds, commitment, point, result_blind)
    }

    #[test]
    fn hyrax_prove_commits_to_claimed_evaluation() {
        let (pcs, evaluations, blinds, commitment, point, result_blind) = example_statement();
        let mut transcript = Transcript::new(b"hyrax-opening");
        let mut rng = test_rng();

        let proof = pcs.prove_with_transcript(
            &commitment,
            &evaluations,
            &blinds,
            &point,
            &result_blind,
            &mut transcript,
            &mut rng,
        );

        let expected_value = mle_from_evaluations(evaluations).eval(&point);
        assert_eq!(
            proof.result_com,
            pcs.commit_scalar(&expected_value, &result_blind)
        );
    }

    #[test]
    fn hyrax_verify_accepts_valid_opening() {
        let (pcs, evaluations, blinds, commitment, point, result_blind) = example_statement();
        let mut prover_transcript = Transcript::new(b"hyrax-opening");
        let mut rng = test_rng();
        let proof = pcs.prove_with_transcript(
            &commitment,
            &evaluations,
            &blinds,
            &point,
            &result_blind,
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"hyrax-opening");
        assert!(pcs.verify_with_transcript(&commitment, &point, &proof, &mut verifier_transcript));
    }

    #[test]
    fn hyrax_verify_rejects_different_point() {
        let (pcs, evaluations, blinds, commitment, point, result_blind) = example_statement();
        let mut prover_transcript = Transcript::new(b"hyrax-opening");
        let mut rng = test_rng();
        let proof = pcs.prove_with_transcript(
            &commitment,
            &evaluations,
            &blinds,
            &point,
            &result_blind,
            &mut prover_transcript,
            &mut rng,
        );
        let other_point = [F::from(2), F::from(3), F::from(5), F::from(8)];

        let mut verifier_transcript = Transcript::new(b"hyrax-opening");
        assert!(!pcs.verify_with_transcript(
            &commitment,
            &other_point,
            &proof,
            &mut verifier_transcript
        ));
    }

    #[test]
    fn hyrax_verify_rejects_tampered_proof() {
        let (pcs, evaluations, blinds, commitment, point, result_blind) = example_statement();
        let mut prover_transcript = Transcript::new(b"hyrax-opening");
        let mut rng = test_rng();
        let mut proof = pcs.prove_with_transcript(
            &commitment,
            &evaluations,
            &blinds,
            &point,
            &result_blind,
            &mut prover_transcript,
            &mut rng,
        );
        proof.schnorr_res.0 += F::from(1);

        let mut verifier_transcript = Transcript::new(b"hyrax-opening");
        assert!(!pcs.verify_with_transcript(&commitment, &point, &proof, &mut verifier_transcript));
    }
}
