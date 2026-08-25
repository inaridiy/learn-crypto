//! Hyrax の square-root multilinear polynomial commitment scheme。
//!
//! 評価表を column-major な正方行列 `T` と見なし、各行を Pedersen commit する。
//! 点 `r = (r_L, r_R)` での opening は、行 commitment を `eq(r_L)` で射影した後、
//! `sigma::InnerProductProof` で `eq(r_R)` との内積関係を証明する。

use ark_ec::{
    CurveGroup,
    hashing::{HashToCurve, HashToCurveError},
};
use ark_std::rand::{CryptoRng, Rng};

use crate::{
    primitive::{
        EqPoly, MultilinearPoly, ScalarPedersen, Transcript, VectorPedersen, column_major_row,
        inner_product,
    },
    protocol::sigma::InnerProductProof,
};

#[derive(Clone, Debug)]
pub struct HyraxPCS<G: CurveGroup> {
    /// `T` の一行を commit する multi-commitment key。
    pub rows: VectorPedersen<G>,
    /// 評価値を commit する scalar commitment key。
    pub scalar: ScalarPedersen<G>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HyraxPCSCommitment<G: CurveGroup> {
    pub rows: Vec<G>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HyraxPCSProof<G: CurveGroup> {
    /// `m(r)` への hiding Pedersen commitment。
    pub result_com: G,
    pub ipa: InnerProductProof<G>,
}

impl<G: CurveGroup> HyraxPCS<G> {
    /// `square_size = 2^(num_vars/2)` 個の行・列を持つ square-root PCS を設定する。
    ///
    /// 論文の commitment key と同じく、vector/scalar commitment は blind generator
    /// `h` を共有し、値側の generators だけを分離する。
    pub fn setup<H: HashToCurve<G>>(
        domain: &[u8],
        square_size: usize,
    ) -> Result<Self, HashToCurveError> {
        assert!(
            square_size.is_power_of_two(),
            "Hyrax square size must be a non-zero power of two"
        );

        let rows = VectorPedersen::setup::<H>(domain, square_size)?;
        let scalar = ScalarPedersen::setup::<H>(domain)?;
        debug_assert_eq!(rows.blind, scalar.blind);

        Ok(Self { rows, scalar })
    }

    pub fn commit<P>(&self, poly: &P, com_blinds: &[G::ScalarField]) -> HyraxPCSCommitment<G>
    where
        P: MultilinearPoly<G::ScalarField> + ?Sized,
    {
        self.assert_poly_shape(poly);
        assert_eq!(
            com_blinds.len(),
            self.rows.len(),
            "one commitment blind is required per matrix row"
        );

        let evaluations = poly.to_evaluations();
        assert_eq!(
            evaluations.len(),
            self.rows.len() * self.rows.len(),
            "polynomial evaluation table does not match the Hyrax matrix"
        );
        let rows = (0..self.rows.len())
            .map(|row| {
                let values =
                    column_major_row(&evaluations, self.rows.len(), row).collect::<Vec<_>>();
                self.rows.commit(&values, &com_blinds[row])
            })
            .collect();

        HyraxPCSCommitment { rows }
    }

    fn row_vars(&self) -> usize {
        debug_assert!(self.rows.len().is_power_of_two());
        self.rows.len().ilog2() as usize
    }

    fn append_statement(
        transcript: &mut Transcript,
        commitment: &HyraxPCSCommitment<G>,
        point: &[G::ScalarField],
        result_com: &G,
    ) {
        transcript.append_usize(b"hyrax-square-size", commitment.rows.len());
        for row in &commitment.rows {
            transcript.append_serializable(b"hyrax-row-commitment", row);
        }
        for coordinate in point {
            transcript.append_serializable(b"hyrax-point", coordinate);
        }
        transcript.append_serializable(b"hyrax-result-commitment", result_com);
    }

    /// 点 `point` での評価証明を作る(PCS の `Open`)。
    #[allow(clippy::too_many_arguments)]
    pub fn open<P>(
        &self,
        commitment: &HyraxPCSCommitment<G>,
        poly: &P,
        com_blinds: &[G::ScalarField],
        point: &[G::ScalarField],
        result_blind: &G::ScalarField,
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> HyraxPCSProof<G>
    where
        P: MultilinearPoly<G::ScalarField> + Clone,
    {
        self.assert_opening_shape(commitment, poly, com_blinds, point);

        let split = self.row_vars();
        let mut row_reduced_poly = poly.clone();
        for &coordinate in &point[..split] {
            row_reduced_poly.fold(coordinate);
        }
        let row_reduced_evaluations = row_reduced_poly.to_evaluations();
        assert_eq!(
            row_reduced_evaluations.len(),
            self.rows.len(),
            "partially evaluated polynomial does not match a Hyrax row"
        );
        let (row_weights, column_weights) = (
            EqPoly::new(point[..split].to_vec()).table(),
            EqPoly::new(point[split..].to_vec()).table(),
        );
        let row_reduced_blind = inner_product(com_blinds, &row_weights);
        let row_reduced_commitment = inner_product(&commitment.rows, &row_weights);
        let result = inner_product(&row_reduced_evaluations, &column_weights);
        debug_assert_eq!(result, poly.eval(point));
        let result_com = self.scalar.commit(&result, result_blind);

        Self::append_statement(transcript, commitment, point, &result_com);
        let ipa = InnerProductProof::prove(
            &self.rows,
            &self.scalar,
            &row_reduced_commitment,
            &result_com,
            &row_reduced_evaluations,
            &row_reduced_blind,
            &column_weights,
            result_blind,
            transcript,
            rng,
        );

        HyraxPCSProof { result_com, ipa }
    }

    pub fn verify(
        &self,
        commitment: &HyraxPCSCommitment<G>,
        point: &[G::ScalarField],
        proof: &HyraxPCSProof<G>,
        transcript: &mut Transcript,
    ) -> bool {
        let split = self.row_vars();
        if commitment.rows.len() != self.rows.len()
            || point.len() != 2 * split
            || self.rows.blind != self.scalar.blind
        {
            return false;
        }

        let (row_weights, column_weights) = (
            EqPoly::new(point[..split].to_vec()).table(),
            EqPoly::new(point[split..].to_vec()).table(),
        );
        let row_reduced_commitment = inner_product(&commitment.rows, &row_weights);

        Self::append_statement(transcript, commitment, point, &proof.result_com);
        proof.ipa.verify(
            &self.rows,
            &self.scalar,
            &column_weights,
            &row_reduced_commitment,
            &proof.result_com,
            transcript,
        )
    }

    fn assert_poly_shape<P>(&self, poly: &P)
    where
        P: MultilinearPoly<G::ScalarField> + ?Sized,
    {
        assert_eq!(
            poly.vars(),
            2 * self.row_vars(),
            "square Hyrax requires exactly 2*log2(square_size) variables"
        );
    }

    fn assert_opening_shape<P>(
        &self,
        commitment: &HyraxPCSCommitment<G>,
        poly: &P,
        com_blinds: &[G::ScalarField],
        point: &[G::ScalarField],
    ) where
        P: MultilinearPoly<G::ScalarField> + ?Sized,
    {
        self.assert_poly_shape(poly);
        assert_eq!(
            commitment.rows.len(),
            self.rows.len(),
            "commitment row count does not match the PCS"
        );
        assert_eq!(
            com_blinds.len(),
            self.rows.len(),
            "one commitment blind is required per matrix row"
        );
        assert_eq!(
            point.len(),
            poly.vars(),
            "evaluation point dimension does not match the polynomial"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::{HyraxPCS, HyraxPCSCommitment, HyraxPCSProof};
    use crate::primitive::{DenseMultilinearPoly, MultilinearPoly, Transcript};
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

    #[allow(clippy::type_complexity)]
    fn example() -> (
        HyraxPCS<G1Projective>,
        DenseMultilinearPoly<F>,
        Vec<F>,
        HyraxPCSCommitment<G1Projective>,
        Vec<F>,
        F,
    ) {
        let pcs = HyraxPCS::setup::<G1Hasher>(b"hyrax-test", 4).unwrap();
        let poly =
            DenseMultilinearPoly::new((1..=16).map(|i| F::from((i * i) as u64)).collect(), 4);
        let blinds = [3, 5, 7, 11].map(F::from).to_vec();
        let commitment = pcs.commit(&poly, &blinds);
        let point = [2, 3, 5, 7].map(F::from).to_vec();
        let result_blind = F::from(13);
        (pcs, poly, blinds, commitment, point, result_blind)
    }

    #[allow(clippy::type_complexity)]
    fn prove_example() -> (
        HyraxPCS<G1Projective>,
        DenseMultilinearPoly<F>,
        HyraxPCSCommitment<G1Projective>,
        Vec<F>,
        F,
        HyraxPCSProof<G1Projective>,
    ) {
        let (pcs, poly, blinds, commitment, point, result_blind) = example();
        let mut transcript = Transcript::new(b"hyrax-opening");
        let proof = pcs.open(
            &commitment,
            &poly,
            &blinds,
            &point,
            &result_blind,
            &mut transcript,
            &mut test_rng(),
        );
        (pcs, poly, commitment, point, result_blind, proof)
    }

    #[test]
    fn commitment_uses_the_papers_column_major_matrix_layout() {
        let (pcs, poly, blinds, commitment, _, _) = example();
        let evals = poly.to_evaluations();
        for row in 0..4 {
            let expected_row = [evals[row], evals[row + 4], evals[row + 8], evals[row + 12]];
            assert_eq!(
                commitment.rows[row],
                pcs.rows.commit(&expected_row, &blinds[row])
            );
        }
    }

    #[test]
    fn proof_commits_to_the_mle_evaluation_and_verifies() {
        let (pcs, poly, commitment, point, result_blind, proof) = prove_example();
        assert_eq!(
            proof.result_com,
            pcs.scalar.commit(&poly.eval(&point), &result_blind)
        );

        let mut transcript = Transcript::new(b"hyrax-opening");
        assert!(pcs.verify(&commitment, &point, &proof, &mut transcript));
    }

    #[test]
    fn dense_and_sparse_polynomials_have_identical_commitments_and_openings() {
        let (pcs, dense, blinds, dense_commitment, point, result_blind) = example();
        let sparse = dense.to_sparse();
        let sparse_commitment = pcs.commit(&sparse, &blinds);
        assert_eq!(dense_commitment, sparse_commitment);

        let mut prover_transcript = Transcript::new(b"hyrax-opening");
        let proof = pcs.open(
            &sparse_commitment,
            &sparse,
            &blinds,
            &point,
            &result_blind,
            &mut prover_transcript,
            &mut test_rng(),
        );
        let mut verifier_transcript = Transcript::new(b"hyrax-opening");
        assert!(pcs.verify(&sparse_commitment, &point, &proof, &mut verifier_transcript));
    }

    #[test]
    fn verifier_rejects_changed_point_commitment_and_proof() {
        let (pcs, _, commitment, mut point, _, proof) = prove_example();
        point[3] += F::from(1);
        let mut transcript = Transcript::new(b"hyrax-opening");
        assert!(!pcs.verify(&commitment, &point, &proof, &mut transcript));

        let (pcs, _, mut commitment, point, _, proof) = prove_example();
        commitment.rows[0] += pcs.rows.generators[0];
        let mut transcript = Transcript::new(b"hyrax-opening");
        assert!(!pcs.verify(&commitment, &point, &proof, &mut transcript));

        let (pcs, _, commitment, point, _, mut proof) = prove_example();
        proof.result_com += pcs.scalar.generator;
        let mut transcript = Transcript::new(b"hyrax-opening");
        assert!(!pcs.verify(&commitment, &point, &proof, &mut transcript));

        let (pcs, _, commitment, point, _, mut proof) = prove_example();
        proof.ipa.dot_product.z_beta += F::from(1);
        let mut transcript = Transcript::new(b"hyrax-opening");
        assert!(!pcs.verify(&commitment, &point, &proof, &mut transcript));
    }

    #[test]
    fn zero_variable_polynomial_uses_the_length_one_ipa_base_case() {
        let pcs = HyraxPCS::setup::<G1Hasher>(b"hyrax-constant", 1).unwrap();
        let poly = DenseMultilinearPoly::new(vec![F::from(42)], 0);
        let blinds = [F::from(7)];
        let commitment = pcs.commit(&poly, &blinds);
        let mut prover_transcript = Transcript::new(b"hyrax-opening");
        let proof = pcs.open(
            &commitment,
            &poly,
            &blinds,
            &[],
            &F::from(11),
            &mut prover_transcript,
            &mut test_rng(),
        );
        assert!(proof.ipa.minus.is_empty());

        let mut verifier_transcript = Transcript::new(b"hyrax-opening");
        assert!(pcs.verify(&commitment, &[], &proof, &mut verifier_transcript));
    }

    #[test]
    fn verifier_rejects_invalid_public_shapes_without_panicking() {
        let (pcs, _, commitment, point, _, mut proof) = prove_example();
        proof.ipa.plus.pop();
        let mut transcript = Transcript::new(b"hyrax-opening");
        assert!(!pcs.verify(&commitment, &point, &proof, &mut transcript));

        let mut transcript = Transcript::new(b"hyrax-opening");
        assert!(!pcs.verify(
            &HyraxPCSCommitment { rows: vec![] },
            &point,
            &proof,
            &mut transcript
        ));
    }
}
