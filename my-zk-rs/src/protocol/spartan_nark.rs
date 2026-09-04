//! Spartan NARK。R1CS instance $A z \circ B z = C z$ の充足性を、
//! witness MLE $\tilde{W}$ への Hyrax commitment と sumcheck で証明し、
//! Fiat--Shamir 変換で非対話化する。zero-knowledge compiler は含まない。
//!
//! `z = [witness..., 0..., 1, io...,  0...]`(Spartan 配置)のうち、
//! prover が commit するのは witness 側の半分だけ。public 側 $\widetilde{(1, io)}$ は
//! verifier が自力で評価できる。

use ark_ec::{
    hashing::{HashToCurve, HashToCurveError},
    CurveGroup,
};
use ark_ff::{UniformRand, Zero};
use ark_std::rand::{CryptoRng, Rng};

use crate::{
    primitive::{DenseMultilinearPoly, EqPoly, Matrix, MultilinearPoly, SpartanR1CS, Transcript},
    protocol::{
        hyrax::{HyraxPCS, HyraxPCSCommitment, HyraxPCSProof},
        sumcheck::{SumcheckOutput, SumcheckProof},
    },
};

pub struct SpartanNarkInstance<G: CurveGroup, M: Matrix<G::ScalarField>> {
    pcs: HyraxPCS<G>,
    r1cs: SpartanR1CS<G::ScalarField, M>,
}

pub struct SpartanNarkProof<G: CurveGroup> {
    /// $\tilde{W}$ への Hyrax commitment(prover の最初の message)。
    pub witness_com: HyraxPCSCommitment<G>,
    /// 1 本目の sumcheck:
    /// $0 = \sum_x \mathrm{eq}(\tau, x) (\tilde{Az}(x) \tilde{Bz}(x) - \tilde{Cz}(x))$ の zero-check。
    pub sc1: SumcheckProof<G::ScalarField>,
    /// Sumcheck の終点 $r_x$ での claims $(v_A, v_B, v_C) = (\tilde{Az}(r_x), \tilde{Bz}(r_x), \tilde{Cz}(r_x))$。
    /// verifier は $\mathrm{eq}(\tau, r_x) (v_A v_B - v_C)$ が sumcheck の最終 claim と
    /// 一致するかを確認する。
    pub va: G::ScalarField,
    pub vb: G::ScalarField,
    pub vc: G::ScalarField,
    /// 2 本目の sumcheck:
    /// $\rho_A v_A + \rho_B v_B + \rho_C v_C$ を
    /// $\sum_y (\rho_A \tilde A(r_x,y) + \rho_B \tilde B(r_x,y)
    /// + \rho_C \tilde C(r_x,y)) \tilde Z(y)$ に還元する。
    pub sc2: SumcheckProof<G::ScalarField>,
    /// 2 本目の終点 $r_y$ の witness 側での $\tilde W$ の評価値。
    pub witness_eval: G::ScalarField,
    /// `witness_eval` が最初の commitment の評価であることを示す Hyrax opening。
    pub witness_opening: HyraxPCSProof<G>,
}

impl<G: CurveGroup, M: Matrix<G::ScalarField>> SpartanNarkInstance<G, M> {
    pub fn encode<H: HashToCurve<G>>(
        r1cs: SpartanR1CS<G::ScalarField, M>,
    ) -> Result<Self, HashToCurveError> {
        // commit するのは $s$ 変数の witness MLE $\tilde{W}$。
        let pcs = HyraxPCS::setup::<H>(b"spartan_nark_pcs", r1cs.half_vars())?;
        Ok(Self { pcs, r1cs })
    }

    /// 行変数を `rx` に固定した行列 MLE の、列側 Boolean hypercube 上での評価表。
    fn row_reduced_matrix(
        matrix: &M,
        rx: &[G::ScalarField],
    ) -> DenseMultilinearPoly<G::ScalarField> {
        DenseMultilinearPoly::from_evals(matrix.vec_mul(&EqPoly::new(rx.to_vec()).to_evals()))
    }

    /// 証明対象の statement(R1CS と io)を transcript に bind する。
    fn append_statement(&self, io: &[G::ScalarField], transcript: &mut Transcript) {
        let structure = &self.r1cs.structure;
        transcript.append_usize(b"spartan-num-constraints", structure.num_constraints);
        transcript.append_usize(b"spartan-num-io", structure.num_io);
        transcript.append_usize(b"spartan-num-witness", structure.num_witness);
        transcript.append_matrix(b"spartan-matrix-a", &self.r1cs.a);
        transcript.append_matrix(b"spartan-matrix-b", &self.r1cs.b);
        transcript.append_matrix(b"spartan-matrix-c", &self.r1cs.c);
        for value in io {
            transcript.append_serializable(b"spartan-io", value);
        }
    }

    pub fn prove(
        &self,
        io: &[G::ScalarField],
        witness: &[G::ScalarField],
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> SpartanNarkProof<G> {
        let z = self.r1cs.assignment(io, witness);
        assert!(self.r1cs.is_sat(&z), "assignment does not satisfy the R1CS");

        self.append_statement(io, transcript);

        // Prover の最初の message: 行ごとに blind を引いた $\tilde{W}$ への
        // hiding commitment。以降の challenge に bind するため transcript に積む。
        let witness_mle = self.r1cs.witness_mle(witness);
        let com_blinds = (0..self.pcs.num_rows())
            .map(|_| G::ScalarField::rand(rng))
            .collect::<Vec<_>>();
        let witness_com = self.pcs.commit(&witness_mle, &com_blinds);
        transcript.append_serializable(b"spartan-witness-commitment", &witness_com);

        // Step 1: verifier の challenge $\tau$ を導出し、
        // $\sum_x \mathrm{eq}(\tau, x) (\tilde{Az}(x) \tilde{Bz}(x) - \tilde{Cz}(x)) = 0$ を sumcheck する。
        let constraint_vars = self.r1cs.constraint_vars();
        let tau = (0..constraint_vars)
            .map(|_| transcript.challenge_field::<G::ScalarField>(b"spartan-tau"))
            .collect::<Vec<_>>();
        let eq = EqPoly::new(tau).to_dense_mlp();
        let az = DenseMultilinearPoly::from_evals(self.r1cs.a.mul_vec(&z));
        let bz = DenseMultilinearPoly::from_evals(self.r1cs.b.mul_vec(&z));
        let cz = DenseMultilinearPoly::from_evals(self.r1cs.c.mul_vec(&z));
        let SumcheckOutput {
            proof: sc1,
            point: rx,
            final_evals,
        } = SumcheckProof::prove(
            vec![eq, az, bz, cz],
            G::ScalarField::zero(),
            3,
            |v| v[0] * (v[1] * v[2] - v[3]),
            transcript,
        );

        let (va, vb, vc) = (final_evals[1], final_evals[2], final_evals[3]);
        transcript.append_serializable(b"spartan-claim-a", &va);
        transcript.append_serializable(b"spartan-claim-b", &vb);
        transcript.append_serializable(b"spartan-claim-c", &vc);

        // Step 2: 3 claims を random linear combination し、行列の行変数を
        // r_x に固定した上で、列変数 y に関する一つの sumcheck にまとめる。
        let rho_a = transcript.challenge_field::<G::ScalarField>(b"spartan-rho-a");
        let rho_b = transcript.challenge_field::<G::ScalarField>(b"spartan-rho-b");
        let rho_c = transcript.challenge_field::<G::ScalarField>(b"spartan-rho-c");
        let sc2_claim = rho_a * va + rho_b * vb + rho_c * vc;
        let a_rx = Self::row_reduced_matrix(&self.r1cs.a, &rx);
        let b_rx = Self::row_reduced_matrix(&self.r1cs.b, &rx);
        let c_rx = Self::row_reduced_matrix(&self.r1cs.c, &rx);
        let z_mle = DenseMultilinearPoly::from_evals(z);
        let SumcheckOutput {
            proof: sc2,
            point: ry,
            ..
        } = SumcheckProof::prove(
            vec![a_rx, b_rx, c_rx, z_mle],
            sc2_claim,
            2,
            |v| (rho_a * v[0] + rho_b * v[1] + rho_c * v[2]) * v[3],
            transcript,
        );

        // Step 2 の終点 $r_y$ の witness 側 $s$ 変数で $\tilde{W}$ を開く。
        let witness_point = &ry[..self.r1cs.half_vars()];
        let witness_eval = witness_mle.eval(witness_point);
        transcript.append_serializable(b"spartan-witness-evaluation", &witness_eval);
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

        SpartanNarkProof {
            witness_com,
            sc1,
            va,
            vb,
            vc,
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
        proof: &SpartanNarkProof<G>,
        transcript: &mut Transcript,
    ) -> bool {
        if io.len() != self.r1cs.structure.num_io {
            return false;
        }

        self.append_statement(io, transcript);
        transcript.append_serializable(b"spartan-witness-commitment", &proof.witness_com);

        let constraint_vars = self.r1cs.constraint_vars();
        let tau = (0..constraint_vars)
            .map(|_| transcript.challenge_field::<G::ScalarField>(b"spartan-tau"))
            .collect::<Vec<_>>();
        let Some((sc1_final, rx)) = proof.sc1.verify(G::ScalarField::zero(), 3, transcript) else {
            return false;
        };
        if rx.len() != constraint_vars
            || sc1_final != EqPoly::new(tau).eval(&rx) * (proof.va * proof.vb - proof.vc)
        {
            return false;
        }

        transcript.append_serializable(b"spartan-claim-a", &proof.va);
        transcript.append_serializable(b"spartan-claim-b", &proof.vb);
        transcript.append_serializable(b"spartan-claim-c", &proof.vc);
        let rho_a = transcript.challenge_field::<G::ScalarField>(b"spartan-rho-a");
        let rho_b = transcript.challenge_field::<G::ScalarField>(b"spartan-rho-b");
        let rho_c = transcript.challenge_field::<G::ScalarField>(b"spartan-rho-c");
        let sc2_claim = rho_a * proof.va + rho_b * proof.vb + rho_c * proof.vc;
        let Some((sc2_final, ry)) = proof.sc2.verify(sc2_claim, 2, transcript) else {
            return false;
        };
        if ry.len() != self.r1cs.vars() {
            return false;
        }

        let witness_point = &ry[..self.r1cs.half_vars()];
        transcript.append_serializable(b"spartan-witness-evaluation", &proof.witness_eval);
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

        let z_eval = self.r1cs.assignment_eval(io, &ry, proof.witness_eval);
        let matrix_eval = rho_a * self.r1cs.a.eval_mle(&rx, &ry)
            + rho_b * self.r1cs.b.eval_mle(&rx, &ry)
            + rho_c * self.r1cs.c.eval_mle(&rx, &ry);

        sc2_final == matrix_eval * z_eval
    }
}

#[cfg(test)]
mod tests {
    use super::SpartanNarkInstance;
    use crate::primitive::{
        DenseMatrix, DenseMultilinearPoly, R1CSStructure, SparseMatrix, SpartanR1CS, Transcript,
        R1CS,
    };
    use ark_bls12_381::{g1::Config, Fr as F, G1Projective};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::{field_hashers::DefaultFieldHasher, UniformRand, Zero};
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    /// x * y = t, t * y = out。io = [x, out], witness = [y, t] で
    /// dense 配置は z = [1, x, out, y, t]。制約 2 行なので sumcheck が 1 round 回る。
    fn square_mul_instance() -> SpartanNarkInstance<G1Projective, SparseMatrix<F>> {
        let r1cs = R1CS::<F, _>::new(
            DenseMatrix::from_usize([[0, 1, 0, 0, 0], [0, 0, 0, 0, 1]]),
            DenseMatrix::from_usize([[0, 0, 0, 1, 0], [0, 0, 0, 1, 0]]),
            DenseMatrix::from_usize([[0, 0, 0, 0, 1], [0, 0, 1, 0, 0]]),
            R1CSStructure::new(2, 2, 2),
        );
        SpartanNarkInstance::encode::<G1Hasher>(SpartanR1CS::from(&r1cs)).unwrap()
    }

    fn example_io_witness() -> ([F; 2], [F; 2]) {
        // x = 3, y = 5, t = 15, out = 75
        ([F::from(3), F::from(75)], [F::from(5), F::from(15)])
    }

    #[test]
    fn prove_commits_to_the_witness_mle() {
        let instance = square_mul_instance();
        let (io, witness) = example_io_witness();

        let mut transcript = Transcript::new(b"spartan-nark-test");
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
    fn complete_proof_verifies_and_replays_the_transcript() {
        let instance = square_mul_instance();
        let (io, witness) = example_io_witness();

        let mut prover_transcript = Transcript::new(b"spartan-nark-test");
        let proof = instance.prove(
            &io,
            &witness,
            &mut prover_transcript,
            &mut StdRng::seed_from_u64(42),
        );

        // Verifier は statement から transcript 全体を replay する。
        let mut transcript = Transcript::new(b"spartan-nark-test");
        assert!(instance.verify(&io, &proof, &mut transcript));
        assert_eq!(
            transcript.challenge_field::<F>(b"probe"),
            prover_transcript.challenge_field::<F>(b"probe")
        );
    }

    #[test]
    fn verifier_rejects_changed_statement_and_proof_values() {
        let instance = square_mul_instance();
        let (io, witness) = example_io_witness();

        let prove = || {
            let mut transcript = Transcript::new(b"spartan-nark-test");
            instance.prove(
                &io,
                &witness,
                &mut transcript,
                &mut StdRng::seed_from_u64(42),
            )
        };

        let proof = prove();
        let mut transcript = Transcript::new(b"spartan-nark-test");
        assert!(!instance.verify(&[F::from(4), F::from(75)], &proof, &mut transcript));

        let mut proof = prove();
        proof.va += F::from(1);
        let mut transcript = Transcript::new(b"spartan-nark-test");
        assert!(!instance.verify(&io, &proof, &mut transcript));

        let mut proof = prove();
        proof.witness_eval += F::from(1);
        let mut transcript = Transcript::new(b"spartan-nark-test");
        assert!(!instance.verify(&io, &proof, &mut transcript));
    }

    #[test]
    fn odd_witness_dimension_uses_a_rectangular_hyrax_matrix() {
        // witness * 1 = out。half_len = 2, half_vars = 1 なので、
        // Hyrax の行列は 1 行 2 列になる。
        let r1cs = R1CS::<F, _>::new(
            DenseMatrix::from_usize([[0, 0, 1]]),
            DenseMatrix::from_usize([[1, 0, 0]]),
            DenseMatrix::from_usize([[0, 1, 0]]),
            R1CSStructure::new(1, 1, 1),
        );
        let instance =
            SpartanNarkInstance::<G1Projective, _>::encode::<G1Hasher>(SpartanR1CS::from(&r1cs))
                .unwrap();
        assert_eq!(instance.r1cs.half_vars(), 1);
        assert_eq!(instance.pcs.num_vars, 1);
        assert_eq!(instance.pcs.num_rows(), 1);
        assert_eq!(instance.pcs.rows.len(), 2);

        let io = [F::from(9)];
        let witness = [F::from(9)];
        let mut prover_transcript = Transcript::new(b"spartan-nark-odd-test");
        let proof = instance.prove(
            &io,
            &witness,
            &mut prover_transcript,
            &mut StdRng::seed_from_u64(7),
        );
        let mut verifier_transcript = Transcript::new(b"spartan-nark-odd-test");
        assert!(instance.verify(&io, &proof, &mut verifier_transcript));
    }

    #[test]
    #[should_panic(expected = "assignment does not satisfy the R1CS")]
    fn prove_rejects_unsatisfying_assignments() {
        let instance = square_mul_instance();
        let mut transcript = Transcript::new(b"spartan-nark-test");
        instance.prove(
            &[F::from(3), F::from(74)],
            &[F::from(5), F::from(15)],
            &mut transcript,
            &mut StdRng::seed_from_u64(42),
        );
    }
}
