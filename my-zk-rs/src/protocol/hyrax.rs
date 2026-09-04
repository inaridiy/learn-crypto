//! Hyrax の square-root multilinear polynomial commitment scheme。
//!
//! $n$ 変数の評価表を column-major な $2^{\lfloor n/2 \rfloor} \times 2^{\lceil n/2 \rceil}$
//! 行列 `T` と見なし、各行を Pedersen commit する。$n$ が奇数なら列が行の 2 倍になる。
//! 点 `r = (r_L, r_R)` での opening は、行 commitment を `eq(r_L)` で射影した後、
//! `sigma::InnerProductProof` で `eq(r_R)` との内積関係を証明する。

use crate::{
    primitive::{
        ColumnMajorMatrix, EqPoly, Matrix, MultilinearPoly, ScalarPedersen, Transcript,
        VectorPedersen, inner_product,
    },
    protocol::sigma::InnerProductProof,
};
use ark_ec::{
    CurveGroup,
    hashing::{HashToCurve, HashToCurveError},
};
use ark_serialize::CanonicalSerialize;
use ark_std::rand::{CryptoRng, Rng};

#[derive(Clone, Debug)]
pub struct HyraxPCS<G: CurveGroup> {
    /// commit できる多項式の変数の数 $n$。
    pub num_vars: usize,
    /// `T` の一行(長さ $2^{\lceil n/2 \rceil}$)を commit する multi-commitment key。
    pub rows: VectorPedersen<G>,
    /// 評価値を commit する scalar commitment key。
    pub scalar: ScalarPedersen<G>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
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
    /// `num_vars` 変数の多項式向けの PCS を設定する。
    ///
    /// 論文の commitment key と同じく、vector/scalar commitment は blind generator
    /// `h` を共有し、値側の generators だけを分離する。
    pub fn setup<H: HashToCurve<G>>(
        domain: &[u8],
        num_vars: usize,
    ) -> Result<Self, HashToCurveError> {
        let col_vars = num_vars - num_vars / 2;
        let rows = VectorPedersen::setup::<H>(domain, 1 << col_vars)?;
        let scalar = ScalarPedersen::setup::<H>(domain)?;
        debug_assert_eq!(rows.blind, scalar.blind);

        Ok(Self {
            num_vars,
            rows,
            scalar,
        })
    }

    /// 行を選ぶ変数の数 $\lfloor n/2 \rfloor$。下位の変数が行、上位の変数が列を選ぶ。
    fn row_vars(&self) -> usize {
        self.num_vars / 2
    }

    /// `T` の行数 $2^{\lfloor n/2 \rfloor}$。commitment と blind の個数。
    pub fn num_rows(&self) -> usize {
        1 << self.row_vars()
    }

    /// 評価表を column-major な行列 `T` と見なす。
    /// $\tilde{f}(r_L, r_R) = \mathrm{eq}(r_L)^\top T \, \mathrm{eq}(r_R)$。
    fn evaluation_matrix<P>(&self, poly: &P) -> ColumnMajorMatrix<G::ScalarField>
    where
        P: MultilinearPoly<G::ScalarField> + ?Sized,
    {
        ColumnMajorMatrix::new(poly.to_dense().into_evals(), self.num_rows())
    }

    pub fn commit<P>(&self, poly: &P, com_blinds: &[G::ScalarField]) -> HyraxPCSCommitment<G>
    where
        P: MultilinearPoly<G::ScalarField> + ?Sized,
    {
        self.assert_poly_shape(poly);
        assert_eq!(
            com_blinds.len(),
            self.num_rows(),
            "one commitment blind is required per matrix row"
        );

        let t = self.evaluation_matrix(poly);
        let rows = (0..t.rows())
            .map(|row| {
                self.rows
                    .commit(&t.row(row).collect::<Vec<_>>(), &com_blinds[row])
            })
            .collect();

        HyraxPCSCommitment { rows }
    }

    fn append_statement(
        transcript: &mut Transcript,
        commitment: &HyraxPCSCommitment<G>,
        point: &[G::ScalarField],
        result_com: &G,
    ) {
        transcript.append_usize(b"hyrax-num-rows", commitment.rows.len());
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
        P: MultilinearPoly<G::ScalarField> + ?Sized,
    {
        self.assert_opening_shape(commitment, poly, com_blinds, point);

        let split = self.row_vars();
        let (row_weights, column_weights) = (
            EqPoly::new(point[..split].to_vec()).to_evals(),
            EqPoly::new(point[split..].to_vec()).to_evals(),
        );
        // 行を eq(r_L) で畳んだ eq(r_L)^T T。その後 eq(r_R) との内積が評価値になる。
        let row_reduced_evaluations = self.evaluation_matrix(poly).vec_mul(&row_weights);
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
        if commitment.rows.len() != self.num_rows()
            || point.len() != self.num_vars
            || self.rows.blind != self.scalar.blind
        {
            return false;
        }

        let (row_weights, column_weights) = (
            EqPoly::new(point[..split].to_vec()).to_evals(),
            EqPoly::new(point[split..].to_vec()).to_evals(),
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
            self.num_vars,
            "polynomial variable count does not match the PCS"
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
            self.num_rows(),
            "commitment row count does not match the PCS"
        );
        assert_eq!(
            com_blinds.len(),
            self.num_rows(),
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
        let evals = poly.evals();
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
    fn odd_variable_count_uses_a_rectangular_matrix() {
        // 3 変数: 行 2^1 = 2、列 2^2 = 4。
        let pcs = HyraxPCS::setup::<G1Hasher>(b"hyrax-odd", 3).unwrap();
        assert_eq!(pcs.num_rows(), 2);
        assert_eq!(pcs.rows.len(), 4);

        let poly = DenseMultilinearPoly::new((1..=8).map(F::from).collect(), 3);
        let blinds = [3, 5].map(F::from).to_vec();
        let commitment = pcs.commit(&poly, &blinds);
        assert_eq!(commitment.rows.len(), 2);

        let point = [2, 3, 5].map(F::from).to_vec();
        let result_blind = F::from(13);
        let mut prover_transcript = Transcript::new(b"hyrax-opening");
        let proof = pcs.open(
            &commitment,
            &poly,
            &blinds,
            &point,
            &result_blind,
            &mut prover_transcript,
            &mut test_rng(),
        );
        assert_eq!(
            proof.result_com,
            pcs.scalar.commit(&poly.eval(&point), &result_blind)
        );

        let mut verifier_transcript = Transcript::new(b"hyrax-opening");
        assert!(pcs.verify(&commitment, &point, &proof, &mut verifier_transcript));
    }

    #[test]
    fn zero_variable_polynomial_uses_the_length_one_ipa_base_case() {
        let pcs = HyraxPCS::setup::<G1Hasher>(b"hyrax-constant", 0).unwrap();
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
