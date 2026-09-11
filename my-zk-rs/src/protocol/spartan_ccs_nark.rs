//! CCS 向けの Spartan NARK。[`SpartanNarkInstance`](super::spartan_nark::SpartanNarkInstance)
//! を R1CS から一般の CCS $f(M_1 z, \ldots, M_t z) = 0$ に拡張したもの。
//! zero-knowledge compiler は含まず、sumcheck の終点での claims は平文で送る。
//!
//! R1CS 版からの変更点:
//!
//! - 1 本目の sumcheck は
//!   $0 = \sum_x \mathrm{eq}(\tau, x) f(\widetilde{M_1 z}(x), \ldots, \widetilde{M_t z}(x))$
//!   の zero-check で、ラウンド多項式の次数は $1 + \deg f$ になる。
//! - 終点 $r_x$ での claims $v_j = \widetilde{M_j z}(r_x)$ は $t$ 個あり、verifier は
//!   $\mathrm{eq}(\tau, r_x) f(v_1, \ldots, v_t)$ が sumcheck の最終 claim と一致するかを確認する。
//! - 2 本目の sumcheck の初期 claim は $\sum_j \rho_j v_j$ で、
//!   $\sum_y \left(\sum_j \rho_j \tilde{M_j}(r_x, y)\right) \tilde{Z}(y)$ に還元する。
//!   これは R1CS 版と同じ形で、行列が 3 個から $t$ 個になっただけ。

use ark_ec::{
    CurveGroup,
    hashing::{HashToCurve, HashToCurveError},
};
use ark_ff::{UniformRand, Zero};
use ark_std::rand::{CryptoRng, Rng};

use crate::{
    primitive::{DenseMultilinearPoly, EqPoly, Matrix, MultilinearPoly, SpartanCCS, Transcript},
    protocol::{
        hyrax::{HyraxPCS, HyraxPCSCommitment, HyraxPCSProof},
        sumcheck::{SumcheckOutput, SumcheckProof},
    },
};

pub struct SpartanCcsNarkInstance<G: CurveGroup, M: Matrix<G::ScalarField>> {
    pub pcs: HyraxPCS<G>,
    pub ccs: SpartanCCS<G::ScalarField, M>,
}

pub struct SpartanCcsNarkProof<G: CurveGroup> {
    /// $\tilde{W}$ への Hyrax commitment(prover の最初の message)。
    pub witness_com: HyraxPCSCommitment<G>,
    /// 1 本目の sumcheck:
    /// $0 = \sum_x \mathrm{eq}(\tau, x) f(\widetilde{M_1 z}(x), \ldots, \widetilde{M_t z}(x))$ の zero-check。
    pub sc1: SumcheckProof<G::ScalarField>,
    /// 終点 $r_x$ での claims $v_j = \widetilde{M_j z}(r_x)$($t$ 個)。
    pub claims: Vec<G::ScalarField>,
    /// 2 本目の sumcheck:
    /// $\sum_j \rho_j v_j$ を $\sum_y (\sum_j \rho_j \tilde{M_j}(r_x, y)) \tilde{Z}(y)$ に還元する。
    pub sc2: SumcheckProof<G::ScalarField>,
    /// 2 本目の終点 $r_y$ の witness 側での $\tilde W$ の評価値。
    pub witness_eval: G::ScalarField,
    /// `witness_eval` が最初の commitment の評価であることを示す Hyrax opening。
    pub witness_opening: HyraxPCSProof<G>,
}

impl<G: CurveGroup, M: Matrix<G::ScalarField>> SpartanCcsNarkInstance<G, M> {
    pub fn encode<H: HashToCurve<G>>(
        ccs: SpartanCCS<G::ScalarField, M>,
    ) -> Result<Self, HashToCurveError> {
        // commit するのは $s$ 変数の witness MLE $\tilde{W}$。
        let pcs = HyraxPCS::setup::<H>(b"spartan_ccs_nark_pcs", ccs.half_vars())?;
        Ok(Self { pcs, ccs })
    }

    /// 1 本目の sumcheck のラウンド多項式の次数 $1 + \deg f$。
    fn sc1_degree(&self) -> usize {
        1 + self.ccs.max_degree()
    }

    /// 行変数を `rx` に固定した行列 MLE の、列側 Boolean hypercube 上での評価表。
    fn row_reduced_matrix(
        matrix: &M,
        rx: &[G::ScalarField],
    ) -> DenseMultilinearPoly<G::ScalarField> {
        DenseMultilinearPoly::from_evals(matrix.vec_mul(&EqPoly::new(rx.to_vec()).to_evals()))
    }

    /// 証明対象の statement(CCS と io)を transcript に bind する。
    fn append_statement(&self, io: &[G::ScalarField], transcript: &mut Transcript) {
        let structure = &self.ccs.structure;
        transcript.append_usize(b"spartan-ccs-num-constraints", structure.num_constraints);
        transcript.append_usize(b"spartan-ccs-num-io", structure.num_io);
        transcript.append_usize(b"spartan-ccs-num-witness", structure.num_witness);
        transcript.append_usize(b"spartan-ccs-num-matrices", self.ccs.matrices.len());
        for matrix in &self.ccs.matrices {
            transcript.append_matrix(b"spartan-ccs-matrix", matrix);
        }
        self.ccs
            .poly
            .append_to_transcript(b"spartan-ccs-constraint-poly", transcript);
        for value in io {
            transcript.append_serializable(b"spartan-ccs-io", value);
        }
    }

    fn append_claims(claims: &[G::ScalarField], transcript: &mut Transcript) {
        for claim in claims {
            transcript.append_serializable(b"spartan-ccs-claim", claim);
        }
    }

    fn challenge_rhos(&self, transcript: &mut Transcript) -> Vec<G::ScalarField> {
        (0..self.ccs.matrices.len())
            .map(|_| transcript.challenge_field::<G::ScalarField>(b"spartan-ccs-rho"))
            .collect()
    }

    /// $\sum_j \rho_j v_j$。
    fn weighted_sum(rhos: &[G::ScalarField], values: &[G::ScalarField]) -> G::ScalarField {
        rhos.iter()
            .zip(values)
            .fold(G::ScalarField::zero(), |acc, (&rho, &v)| acc + rho * v)
    }

    pub fn prove(
        &self,
        io: &[G::ScalarField],
        witness: &[G::ScalarField],
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> SpartanCcsNarkProof<G> {
        let z = self.ccs.assignment(io, witness);
        assert!(self.ccs.is_sat(&z), "assignment does not satisfy the CCS");

        let num_matrices = self.ccs.matrices.len();
        self.append_statement(io, transcript);

        // Prover の最初の message: 行ごとに blind を引いた $\tilde{W}$ への
        // hiding commitment。以降の challenge に bind するため transcript に積む。
        let witness_mle = self.ccs.witness_mle(witness);
        let com_blinds = (0..self.pcs.num_rows())
            .map(|_| G::ScalarField::rand(rng))
            .collect::<Vec<_>>();
        let witness_com = self.pcs.commit(&witness_mle, &com_blinds);
        transcript.append_serializable(b"spartan-ccs-witness-commitment", &witness_com);

        // Step 1: verifier の challenge $\tau$ を導出し、
        // $\sum_x \mathrm{eq}(\tau, x) f(\widetilde{M_1 z}(x), \ldots) = 0$ を sumcheck する。
        let constraint_vars = self.ccs.constraint_vars();
        let tau = (0..constraint_vars)
            .map(|_| transcript.challenge_field::<G::ScalarField>(b"spartan-ccs-tau"))
            .collect::<Vec<_>>();
        let mut sc1_polys = Vec::with_capacity(1 + num_matrices);
        sc1_polys.push(EqPoly::new(tau).to_dense_mlp());
        sc1_polys.extend(
            self.ccs
                .matrices
                .iter()
                .map(|m| DenseMultilinearPoly::from_evals(m.mul_vec(&z))),
        );
        let SumcheckOutput {
            proof: sc1,
            point: rx,
            final_evals,
        } = SumcheckProof::prove(
            sc1_polys,
            G::ScalarField::zero(),
            self.sc1_degree(),
            |v| v[0] * self.ccs.poly.eval(&v[1..]),
            transcript,
        );

        let claims = final_evals[1..].to_vec();
        Self::append_claims(&claims, transcript);

        // Step 2: t 個の claims を random linear combination し、行列の行変数を
        // r_x に固定した上で、列変数 y に関する一つの sumcheck にまとめる。
        let rhos = self.challenge_rhos(transcript);
        let sc2_claim = Self::weighted_sum(&rhos, &claims);
        let mut sc2_polys = self
            .ccs
            .matrices
            .iter()
            .map(|m| Self::row_reduced_matrix(m, &rx))
            .collect::<Vec<_>>();
        sc2_polys.push(DenseMultilinearPoly::from_evals(z));
        let SumcheckOutput {
            proof: sc2,
            point: ry,
            ..
        } = SumcheckProof::prove(
            sc2_polys,
            sc2_claim,
            2,
            |v| Self::weighted_sum(&rhos, &v[..num_matrices]) * v[num_matrices],
            transcript,
        );

        // Step 2 の終点 $r_y$ の witness 側 $s$ 変数で $\tilde{W}$ を開く。
        let witness_point = &ry[..self.ccs.half_vars()];
        let witness_eval = witness_mle.eval(witness_point);
        transcript.append_serializable(b"spartan-ccs-witness-evaluation", &witness_eval);
        // 評価値は verifier が使う公開値なので、結果 commitment の blind は 0。
        let witness_opening = self.pcs.open(
            &witness_com,
            &witness_mle,
            &com_blinds,
            witness_point,
            &G::ScalarField::zero(),
            transcript,
            rng,
        );

        SpartanCcsNarkProof {
            witness_com,
            sc1,
            claims,
            sc2,
            witness_eval,
            witness_opening,
        }
    }

    /// Fiat--Shamir transcript を replay し、Spartan の二つの sumcheck と
    /// witness commitment の opening を検証する。
    pub fn verify(
        &self,
        io: &[G::ScalarField],
        proof: &SpartanCcsNarkProof<G>,
        transcript: &mut Transcript,
    ) -> bool {
        if io.len() != self.ccs.structure.num_io || proof.claims.len() != self.ccs.matrices.len() {
            return false;
        }

        self.append_statement(io, transcript);
        transcript.append_serializable(b"spartan-ccs-witness-commitment", &proof.witness_com);

        let constraint_vars = self.ccs.constraint_vars();
        let tau = (0..constraint_vars)
            .map(|_| transcript.challenge_field::<G::ScalarField>(b"spartan-ccs-tau"))
            .collect::<Vec<_>>();
        let Some((sc1_final, rx)) =
            proof
                .sc1
                .verify(G::ScalarField::zero(), self.sc1_degree(), transcript)
        else {
            return false;
        };
        if rx.len() != constraint_vars
            || sc1_final != EqPoly::new(tau).eval(&rx) * self.ccs.poly.eval(&proof.claims)
        {
            return false;
        }

        Self::append_claims(&proof.claims, transcript);
        let rhos = self.challenge_rhos(transcript);
        let sc2_claim = Self::weighted_sum(&rhos, &proof.claims);
        let Some((sc2_final, ry)) = proof.sc2.verify(sc2_claim, 2, transcript) else {
            return false;
        };
        if ry.len() != self.ccs.vars() {
            return false;
        }

        let witness_point = &ry[..self.ccs.half_vars()];
        transcript.append_serializable(b"spartan-ccs-witness-evaluation", &proof.witness_eval);
        if proof.witness_opening.result_com
            != self
                .pcs
                .scalar
                .commit(&proof.witness_eval, &G::ScalarField::zero())
            || !self.pcs.verify(
                &proof.witness_com,
                witness_point,
                &proof.witness_opening,
                transcript,
            )
        {
            return false;
        }

        let z_eval = self.ccs.assignment_eval(io, &ry, proof.witness_eval);
        let matrix_evals = self
            .ccs
            .matrices
            .iter()
            .map(|m| m.eval_mle(&rx, &ry))
            .collect::<Vec<_>>();

        sc2_final == Self::weighted_sum(&rhos, &matrix_evals) * z_eval
    }
}

#[cfg(test)]
mod tests {
    use super::SpartanCcsNarkInstance;
    use crate::primitive::{
        CCS, DenseMatrix, DenseMultilinearPoly, R1CS, R1CSStructure, SparseMatrix, SparseMvPoly,
        SpartanCCS, Transcript,
    };
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::{One, UniformRand, Zero, field_hashers::DefaultFieldHasher};
    use ark_std::rand::{SeedableRng, rngs::StdRng};
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;
    type Instance = SpartanCcsNarkInstance<G1Projective, SparseMatrix<F>>;

    /// x * y = t, t * y = out を R1CS から CCS に変換したもの。
    /// io = [x, out], witness = [y, t] で dense 配置は z = [1, x, out, y, t]。
    fn square_mul_instance() -> Instance {
        let r1cs = R1CS::<F, _>::new(
            DenseMatrix::from_usize([[0, 1, 0, 0, 0], [0, 0, 0, 0, 1]]),
            DenseMatrix::from_usize([[0, 0, 0, 1, 0], [0, 0, 0, 1, 0]]),
            DenseMatrix::from_usize([[0, 0, 0, 0, 1], [0, 0, 1, 0, 0]]),
            R1CSStructure::new(2, 2, 2),
        );
        let ccs = CCS::from(&r1cs);
        SpartanCcsNarkInstance::encode::<G1Hasher>(SpartanCCS::from(&ccs)).unwrap()
    }

    fn square_mul_io_witness() -> ([F; 2], [F; 2]) {
        // x = 3, y = 5, t = 15, out = 75
        ([F::from(3), F::from(75)], [F::from(5), F::from(15)])
    }

    /// 次数 3 の CCS。io = [x, out, s], witness = [y] で z = [1, x, out, s, y]。
    ///
    /// 行 0: x * x * y + x + y - out = 0、行 1: 1 * 1 * y + 1 + y - s = 0。
    /// 制約多項式は $X_0^2 X_1 + X_0 + X_1 - X_2$。
    fn cubic_instance() -> Instance {
        let ccs = CCS::new(
            vec![
                DenseMatrix::<F, 2, 5>::from_usize([[0, 1, 0, 0, 0], [1, 0, 0, 0, 0]]),
                DenseMatrix::from_usize([[0, 0, 0, 0, 1], [0, 0, 0, 0, 1]]),
                DenseMatrix::from_usize([[0, 0, 1, 0, 0], [0, 0, 0, 1, 0]]),
            ],
            SparseMvPoly::new(vec![
                (F::one(), vec![0, 0, 1]),
                (F::one(), vec![0]),
                (F::one(), vec![1]),
                (-F::one(), vec![2]),
            ]),
            R1CSStructure::new(2, 3, 1),
        );
        SpartanCcsNarkInstance::encode::<G1Hasher>(SpartanCCS::from(&ccs)).unwrap()
    }

    fn cubic_io_witness() -> ([F; 3], [F; 1]) {
        // x = 2, y = 3: out = 4 * 3 + 2 + 3 = 17, s = 3 + 1 + 3 = 7
        ([F::from(2), F::from(17), F::from(7)], [F::from(3)])
    }

    #[test]
    fn prove_commits_to_the_witness_mle() {
        let instance = square_mul_instance();
        let (io, witness) = square_mul_io_witness();

        let mut transcript = Transcript::new(b"spartan-ccs-nark-test");
        let proof = instance.prove(
            &io,
            &witness,
            &mut transcript,
            &mut StdRng::seed_from_u64(42),
        );

        // 同じ seed で blind を再現すると、witness 側の半分 [y, t, 0, 0] への
        // Hyrax commitment と一致する。
        let mut rng = StdRng::seed_from_u64(42);
        let blinds = (0..instance.pcs.num_rows())
            .map(|_| F::rand(&mut rng))
            .collect::<Vec<_>>();
        let expected_mle =
            DenseMultilinearPoly::new(vec![F::from(5), F::from(15), F::zero(), F::zero()], 2);
        assert_eq!(
            proof.witness_com,
            instance.pcs.commit(&expected_mle, &blinds)
        );
    }

    #[test]
    fn r1cs_derived_ccs_proof_verifies_and_replays_the_transcript() {
        let instance = square_mul_instance();
        let (io, witness) = square_mul_io_witness();
        assert_eq!(instance.sc1_degree(), 3);

        let mut prover_transcript = Transcript::new(b"spartan-ccs-nark-test");
        let proof = instance.prove(
            &io,
            &witness,
            &mut prover_transcript,
            &mut StdRng::seed_from_u64(42),
        );
        assert_eq!(proof.claims.len(), 3);

        // Verifier は statement から transcript 全体を replay する。
        let mut transcript = Transcript::new(b"spartan-ccs-nark-test");
        assert!(instance.verify(&io, &proof, &mut transcript));
        assert_eq!(
            transcript.challenge_field::<F>(b"probe"),
            prover_transcript.challenge_field::<F>(b"probe")
        );
    }

    #[test]
    fn cubic_ccs_proof_verifies_with_a_degree_four_sumcheck() {
        let instance = cubic_instance();
        let (io, witness) = cubic_io_witness();
        assert_eq!(instance.sc1_degree(), 4);

        let mut prover_transcript = Transcript::new(b"spartan-ccs-nark-test");
        let proof = instance.prove(
            &io,
            &witness,
            &mut prover_transcript,
            &mut StdRng::seed_from_u64(42),
        );

        let mut transcript = Transcript::new(b"spartan-ccs-nark-test");
        assert!(instance.verify(&io, &proof, &mut transcript));
    }

    #[test]
    fn verifier_rejects_changed_statement_and_proof_values() {
        let instance = cubic_instance();
        let (io, witness) = cubic_io_witness();

        let prove = || {
            let mut transcript = Transcript::new(b"spartan-ccs-nark-test");
            instance.prove(
                &io,
                &witness,
                &mut transcript,
                &mut StdRng::seed_from_u64(42),
            )
        };
        let verify = |io: &[F], proof| {
            let mut transcript = Transcript::new(b"spartan-ccs-nark-test");
            instance.verify(io, &proof, &mut transcript)
        };

        assert!(!verify(&[F::from(2), F::from(18), F::from(7)], prove()));

        let mut proof = prove();
        proof.claims[0] += F::from(1);
        assert!(!verify(&io, proof));

        let mut proof = prove();
        proof.claims.pop();
        assert!(!verify(&io, proof));

        let mut proof = prove();
        proof.witness_eval += F::from(1);
        assert!(!verify(&io, proof));
    }

    #[test]
    fn single_constraint_instance_has_a_zero_round_first_sumcheck() {
        // witness * 1 = out。制約 1 行なので constraint_vars = 0 で 1 本目の sumcheck は
        // ラウンドを持たない。
        let r1cs = R1CS::<F, _>::new(
            DenseMatrix::from_usize([[0, 0, 1]]),
            DenseMatrix::from_usize([[1, 0, 0]]),
            DenseMatrix::from_usize([[0, 1, 0]]),
            R1CSStructure::new(1, 1, 1),
        );
        let instance = Instance::encode::<G1Hasher>(SpartanCCS::from(&CCS::from(&r1cs))).unwrap();
        assert_eq!(instance.ccs.constraint_vars(), 0);

        let io = [F::from(9)];
        let witness = [F::from(9)];
        let mut prover_transcript = Transcript::new(b"spartan-ccs-nark-odd-test");
        let proof = instance.prove(
            &io,
            &witness,
            &mut prover_transcript,
            &mut StdRng::seed_from_u64(7),
        );

        let mut verifier_transcript = Transcript::new(b"spartan-ccs-nark-odd-test");
        assert!(instance.verify(&io, &proof, &mut verifier_transcript));
    }

    #[test]
    #[should_panic(expected = "assignment does not satisfy the CCS")]
    fn prove_rejects_unsatisfying_assignments() {
        let instance = cubic_instance();
        let mut transcript = Transcript::new(b"spartan-ccs-nark-test");
        instance.prove(
            &[F::from(2), F::from(18), F::from(7)],
            &[F::from(3)],
            &mut transcript,
            &mut StdRng::seed_from_u64(42),
        );
    }
}
