//! Spartan NIZK。[`SpartanNarkInstance`](super::spartan_nark::SpartanNarkInstance) と
//! 同じ二段の sumcheck で R1CS の充足性を示すが、prover が送る値をすべて Pedersen
//! commitment に置き換え、sigma protocol で関係を示すことで zero knowledge にする。
//!
//! NARK からの変更点:
//!
//! - 両 sumcheck を [`ZkSumcheckProof`] にする。verifier は claim を commitment としてしか持たない。
//! - 1 本目の終点の claims $(v_A, v_B, v_C)$ と積 $v_A v_B$ は commitment で送り、
//!   [`ProductProof`] / [`KnowledgeProof`] で opening を知っていることを示す。
//!   sumcheck の最終 claim $\mathrm{eq}(\tau, r_x)(v_A v_B - v_C)$ との一致は、
//!   commitment の準同型性で組み立てた $\mathrm{eq}(\tau, r_x)(C_{AB} - C_C)$ と
//!   sumcheck が返した commitment の [`EqualityProof`] で示す。
//! - 2 本目の初期 claim $\rho_A v_A + \rho_B v_B + \rho_C v_C$ も commitment の線形結合で作る。
//! - $\tilde{W}(r_y)$ は Hyrax の hiding な評価 commitment として送り、
//!   $\tilde{Z}(r_y)$ の commitment を $(1 - r_s) C_W + r_s \widetilde{(1, io)}(r_y) G$ で組み立て、
//!   行列の評価値(verifier が自力で計算)を掛けたものと 2 本目の最終 claim の一致を
//!   [`EqualityProof`] で示す。
//!
//! 論文の SpartanNIZK と同じく、行列 MLE $\tilde A, \tilde B, \tilde C$ の評価は verifier が
//! 直接計算する(SPARK による sparse polynomial commitment は使わない)。

use ark_ec::{
    CurveGroup,
    hashing::{HashToCurve, HashToCurveError},
};
use ark_ff::{One, UniformRand, Zero};
use ark_std::rand::{CryptoRng, Rng};

use crate::{
    primitive::{DenseMultilinearPoly, EqPoly, Matrix, MultilinearPoly, SpartanR1CS, Transcript},
    protocol::{
        hyrax::{HyraxPCS, HyraxPCSCommitment, HyraxPCSProof},
        sigma::{EqualityProof, KnowledgeProof, ProductProof},
        zk_sumcheck::{ZkSumcheckKey, ZkSumcheckOutput, ZkSumcheckProof},
    },
};

pub struct SpartanNizkInstance<G: CurveGroup, M: Matrix<G::ScalarField>> {
    pcs: HyraxPCS<G>,
    /// 1 本目(次数 3)の sumcheck の key。
    sc1_key: ZkSumcheckKey<G>,
    /// 2 本目(次数 2)の sumcheck の key。
    sc2_key: ZkSumcheckKey<G>,
    r1cs: SpartanR1CS<G::ScalarField, M>,
}

pub struct SpartanNizkProof<G: CurveGroup> {
    /// $\tilde{W}$ への Hyrax commitment(prover の最初の message)。
    pub witness_com: HyraxPCSCommitment<G>,
    /// 1 本目の ZK sumcheck:
    /// $0 = \sum_x \mathrm{eq}(\tau, x) (\tilde{Az}(x) \tilde{Bz}(x) - \tilde{Cz}(x))$ の zero-check。
    pub sc1: ZkSumcheckProof<G>,
    /// 終点 $r_x$ での claims $(v_A, v_B, v_C)$ と積 $v_A v_B$ の commitment。
    pub va_com: G,
    pub vb_com: G,
    pub vc_com: G,
    pub prod_com: G,
    /// `vc_com` の opening の知識証明。
    pub vc_knowledge: KnowledgeProof<G>,
    /// `prod_com` が `va_com` と `vb_com` の積を隠していることの証明。
    pub product: ProductProof<G>,
    /// $\mathrm{eq}(\tau, r_x)(C_{AB} - C_C)$ と 1 本目の最終 claim commitment の一致。
    pub sc1_equality: EqualityProof<G>,
    /// 2 本目の ZK sumcheck:
    /// $\rho_A v_A + \rho_B v_B + \rho_C v_C$ を
    /// $\sum_y (\rho_A \tilde A(r_x,y) + \rho_B \tilde B(r_x,y)
    /// + \rho_C \tilde C(r_x,y)) \tilde Z(y)$ に還元する。
    pub sc2: ZkSumcheckProof<G>,
    /// $\tilde{W}(r_y)$ の hiding commitment(`result_com`)と、それが `witness_com` の
    /// 評価であることの Hyrax opening。
    pub witness_opening: HyraxPCSProof<G>,
    /// 行列評価値 $\times$ $C_{Z(r_y)}$ と 2 本目の最終 claim commitment の一致。
    pub sc2_equality: EqualityProof<G>,
}

impl<G: CurveGroup, M: Matrix<G::ScalarField>> SpartanNizkInstance<G, M> {
    pub fn encode<H: HashToCurve<G>>(
        r1cs: SpartanR1CS<G::ScalarField, M>,
    ) -> Result<Self, HashToCurveError> {
        // commit するのは $s$ 変数の witness MLE $\tilde{W}$。
        let pcs = HyraxPCS::setup::<H>(b"spartan_nizk_pcs", r1cs.half_vars())?;
        // claim の commitment は Hyrax の評価 commitment と準同型に組み合わせるので、
        // scalar key(と blind generator)を PCS と共有する。
        let sc1_key = ZkSumcheckKey::setup::<H>(b"spartan_nizk_sumcheck_1", 3, pcs.scalar.clone())?;
        let sc2_key = ZkSumcheckKey::setup::<H>(b"spartan_nizk_sumcheck_2", 2, pcs.scalar.clone())?;
        Ok(Self {
            pcs,
            sc1_key,
            sc2_key,
            r1cs,
        })
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

    fn append_claim_commitments(proof_claims: [&G; 4], transcript: &mut Transcript) {
        let [va_com, vb_com, vc_com, prod_com] = proof_claims;
        transcript.append_serializable(b"spartan-claim-a-commitment", va_com);
        transcript.append_serializable(b"spartan-claim-b-commitment", vb_com);
        transcript.append_serializable(b"spartan-claim-c-commitment", vc_com);
        transcript.append_serializable(b"spartan-claim-ab-commitment", prod_com);
    }

    fn challenge_rhos(transcript: &mut Transcript) -> [G::ScalarField; 3] {
        [b"spartan-rho-a", b"spartan-rho-b", b"spartan-rho-c"]
            .map(|label| transcript.challenge_field::<G::ScalarField>(label))
    }

    fn assignment_commitment(
        &self,
        io: &[G::ScalarField],
        ry: &[G::ScalarField],
        witness_eval_com: &G,
    ) -> (G, G::ScalarField) {
        let selector = ry[self.r1cs.half_vars()];
        let witness_weight = G::ScalarField::one() - selector;
        let public_eval = self.r1cs.public_mle(io).eval(&ry[..self.r1cs.half_vars()]);
        (
            *witness_eval_com * witness_weight
                + self.pcs.scalar.generator * (selector * public_eval),
            witness_weight,
        )
    }

    pub fn prove(
        &self,
        io: &[G::ScalarField],
        witness: &[G::ScalarField],
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> SpartanNizkProof<G> {
        let z = self.r1cs.assignment(io, witness);
        assert!(self.r1cs.is_sat(&z), "assignment does not satisfy the R1CS");

        let scalar = &self.pcs.scalar;
        self.append_statement(io, transcript);

        // 行ごとに blind を引いた $\tilde{W}$ への hiding commitment。
        let witness_mle = self.r1cs.witness_mle(witness);
        let com_blinds = (0..self.pcs.num_rows())
            .map(|_| G::ScalarField::rand(rng))
            .collect::<Vec<_>>();
        let witness_com = self.pcs.commit(&witness_mle, &com_blinds);
        transcript.append_serializable(b"spartan-witness-commitment", &witness_com);

        //  challenge $\tau$ を導出し、zero-check を ZK sumcheck する。
        // 初期 claim は公開値 0 なので blind も 0。
        let constraint_vars = self.r1cs.constraint_vars();
        let tau = (0..constraint_vars)
            .map(|_| transcript.challenge_field::<G::ScalarField>(b"spartan-tau"))
            .collect::<Vec<_>>();
        let eq = EqPoly::new(tau).to_dense_mlp();
        let az = DenseMultilinearPoly::from_evals(self.r1cs.a.mul_vec(&z));
        let bz = DenseMultilinearPoly::from_evals(self.r1cs.b.mul_vec(&z));
        let cz = DenseMultilinearPoly::from_evals(self.r1cs.c.mul_vec(&z));
        let ZkSumcheckOutput {
            proof: sc1,
            point: rx,
            final_evals,
            final_claim: sc1_final,
            final_blind: sc1_final_blind,
            final_com: sc1_final_com,
        } = ZkSumcheckProof::prove(
            &self.sc1_key,
            vec![eq, az, bz, cz],
            G::ScalarField::zero(),
            G::ScalarField::zero(),
            |v| v[0] * (v[1] * v[2] - v[3]),
            transcript,
            rng,
        );
        let (eq_rx, va, vb, vc) = (
            final_evals[0],
            final_evals[1],
            final_evals[2],
            final_evals[3],
        );

        // 終点の claims は commitment で送り、sigma protocol で関係を示す。
        let [va_blind, vb_blind, vc_blind, prod_blind] =
            core::array::from_fn(|_| G::ScalarField::rand(rng));
        let va_com = scalar.commit(&va, &va_blind);
        let vb_com = scalar.commit(&vb, &vb_blind);
        let vc_com = scalar.commit(&vc, &vc_blind);
        let prod_com = scalar.commit(&(va * vb), &prod_blind);
        Self::append_claim_commitments([&va_com, &vb_com, &vc_com, &prod_com], transcript);
        let vc_knowledge = KnowledgeProof::prove(scalar, &vc_com, &vc, &vc_blind, transcript, rng);
        let product = ProductProof::prove(
            scalar,
            &va_com,
            &vb_com,
            &prod_com,
            &va,
            &va_blind,
            &vb,
            &vb_blind,
            &prod_blind,
            transcript,
            rng,
        );

        // 1 本目の最終検査: eq(tau, r_x) (v_A v_B - v_C) = sc1 の最終 claim。
        // 左辺は commitment の準同型性で組み立てられる。
        let expected1 = eq_rx * (va * vb - vc);
        let expected1_blind = eq_rx * (prod_blind - vc_blind);
        let expected1_com = (prod_com - vc_com) * eq_rx;
        debug_assert_eq!(expected1, sc1_final);
        let sc1_equality = EqualityProof::prove(
            scalar,
            &expected1_com,
            &sc1_final_com,
            &sc1_final,
            &expected1_blind,
            &sc1_final_blind,
            transcript,
            rng,
        );

        // Step 2: 3 claims を random linear combination し、列変数 y についての
        // 一つの ZK sumcheck にまとめる。初期 claim の commitment も線形結合で作れる。
        let [rho_a, rho_b, rho_c] = Self::challenge_rhos(transcript);
        let sc2_claim = rho_a * va + rho_b * vb + rho_c * vc;
        let sc2_blind = rho_a * va_blind + rho_b * vb_blind + rho_c * vc_blind;
        let a_rx = Self::row_reduced_matrix(&self.r1cs.a, &rx);
        let b_rx = Self::row_reduced_matrix(&self.r1cs.b, &rx);
        let c_rx = Self::row_reduced_matrix(&self.r1cs.c, &rx);
        let z_mle = DenseMultilinearPoly::from_evals(z);
        let ZkSumcheckOutput {
            proof: sc2,
            point: ry,
            final_evals,
            final_claim: sc2_final,
            final_blind: sc2_final_blind,
            final_com: sc2_final_com,
        } = ZkSumcheckProof::prove(
            &self.sc2_key,
            vec![a_rx, b_rx, c_rx, z_mle],
            sc2_claim,
            sc2_blind,
            |v| (rho_a * v[0] + rho_b * v[1] + rho_c * v[2]) * v[3],
            transcript,
            rng,
        );
        let matrix_eval = rho_a * final_evals[0] + rho_b * final_evals[1] + rho_c * final_evals[2];
        let z_eval = final_evals[3];

        // $\tilde{W}(r_y)$ を hiding commitment として開く。
        let witness_point = &ry[..self.r1cs.half_vars()];
        let witness_eval = witness_mle.eval(witness_point);
        let witness_eval_blind = G::ScalarField::rand(rng);
        let witness_opening = self.pcs.open(
            &witness_com,
            &witness_mle,
            &com_blinds,
            witness_point,
            &witness_eval_blind,
            transcript,
            rng,
        );

        // 2 本目の最終検査: (rho . (A, B, C)(r_x, r_y)) * Z(r_y) = sc2 の最終 claim。
        let (z_com, witness_weight) =
            self.assignment_commitment(io, &ry, &witness_opening.result_com);
        debug_assert_eq!(z_eval, self.r1cs.assignment_eval(io, &ry, witness_eval));
        let expected2 = matrix_eval * z_eval;
        let expected2_blind = matrix_eval * witness_weight * witness_eval_blind;
        let expected2_com = z_com * matrix_eval;
        debug_assert_eq!(expected2, sc2_final);
        let sc2_equality = EqualityProof::prove(
            scalar,
            &expected2_com,
            &sc2_final_com,
            &sc2_final,
            &expected2_blind,
            &sc2_final_blind,
            transcript,
            rng,
        );

        SpartanNizkProof {
            witness_com,
            sc1,
            va_com,
            vb_com,
            vc_com,
            prod_com,
            vc_knowledge,
            product,
            sc1_equality,
            sc2,
            witness_opening,
            sc2_equality,
        }
    }

    pub fn verify(
        &self,
        io: &[G::ScalarField],
        proof: &SpartanNizkProof<G>,
        transcript: &mut Transcript,
    ) -> bool {
        if io.len() != self.r1cs.structure.num_io {
            return false;
        }

        let scalar = &self.pcs.scalar;
        self.append_statement(io, transcript);
        transcript.append_serializable(b"spartan-witness-commitment", &proof.witness_com);

        // Step 1: 初期 claim は 0 の commitment(単位元)。
        let constraint_vars = self.r1cs.constraint_vars();
        let tau = (0..constraint_vars)
            .map(|_| transcript.challenge_field::<G::ScalarField>(b"spartan-tau"))
            .collect::<Vec<_>>();
        let Some((sc1_final_com, rx)) = proof.sc1.verify(&self.sc1_key, &G::zero(), transcript)
        else {
            return false;
        };
        if rx.len() != constraint_vars {
            return false;
        }

        Self::append_claim_commitments(
            [&proof.va_com, &proof.vb_com, &proof.vc_com, &proof.prod_com],
            transcript,
        );
        if !proof.vc_knowledge.verify(scalar, &proof.vc_com, transcript)
            || !proof.product.verify(
                scalar,
                &proof.va_com,
                &proof.vb_com,
                &proof.prod_com,
                transcript,
            )
        {
            return false;
        }

        let eq_rx = EqPoly::new(tau).eval(&rx);
        let expected1_com = (proof.prod_com - proof.vc_com) * eq_rx;
        if !proof
            .sc1_equality
            .verify(scalar, &expected1_com, &sc1_final_com, transcript)
        {
            return false;
        }

        // Step 2
        let [rho_a, rho_b, rho_c] = Self::challenge_rhos(transcript);
        let sc2_claim_com = proof.va_com * rho_a + proof.vb_com * rho_b + proof.vc_com * rho_c;
        let Some((sc2_final_com, ry)) = proof.sc2.verify(&self.sc2_key, &sc2_claim_com, transcript)
        else {
            return false;
        };
        if ry.len() != self.r1cs.vars() {
            return false;
        }

        let witness_point = &ry[..self.r1cs.half_vars()];
        if !self.pcs.verify(
            &proof.witness_com,
            witness_point,
            &proof.witness_opening,
            transcript,
        ) {
            return false;
        }

        let (z_com, _) = self.assignment_commitment(io, &ry, &proof.witness_opening.result_com);
        let matrix_eval = rho_a * self.r1cs.a.eval_mle(&rx, &ry)
            + rho_b * self.r1cs.b.eval_mle(&rx, &ry)
            + rho_c * self.r1cs.c.eval_mle(&rx, &ry);
        let expected2_com = z_com * matrix_eval;

        proof
            .sc2_equality
            .verify(scalar, &expected2_com, &sc2_final_com, transcript)
    }
}

#[cfg(test)]
mod tests {
    use super::SpartanNizkInstance;
    use crate::primitive::{
        DenseMatrix, DenseMultilinearPoly, R1CS, R1CSStructure, SparseMatrix, SpartanR1CS,
        Transcript,
    };
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::{UniformRand, Zero, field_hashers::DefaultFieldHasher};
    use ark_std::rand::{SeedableRng, rngs::StdRng};
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    /// x * y = t, t * y = out。io = [x, out], witness = [y, t] で
    /// dense 配置は z = [1, x, out, y, t]。制約 2 行なので sumcheck が 1 round 回る。
    fn square_mul_instance() -> SpartanNizkInstance<G1Projective, SparseMatrix<F>> {
        let r1cs = R1CS::<F, _>::new(
            DenseMatrix::from_usize([[0, 1, 0, 0, 0], [0, 0, 0, 0, 1]]),
            DenseMatrix::from_usize([[0, 0, 0, 1, 0], [0, 0, 0, 1, 0]]),
            DenseMatrix::from_usize([[0, 0, 0, 0, 1], [0, 0, 1, 0, 0]]),
            R1CSStructure::new(2, 2, 2),
        );
        SpartanNizkInstance::encode::<G1Hasher>(SpartanR1CS::from(&r1cs)).unwrap()
    }

    fn example_io_witness() -> ([F; 2], [F; 2]) {
        // x = 3, y = 5, t = 15, out = 75
        ([F::from(3), F::from(75)], [F::from(5), F::from(15)])
    }

    #[test]
    fn prove_commits_to_the_witness_mle() {
        let instance = square_mul_instance();
        let (io, witness) = example_io_witness();

        let mut transcript = Transcript::new(b"spartan-nizk-test");
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

        let mut prover_transcript = Transcript::new(b"spartan-nizk-test");
        let proof = instance.prove(
            &io,
            &witness,
            &mut prover_transcript,
            &mut StdRng::seed_from_u64(42),
        );
        assert_eq!(proof.sc1.poly_coms.len(), 1);
        assert_eq!(proof.sc2.poly_coms.len(), 3);

        // Verifier は statement から transcript 全体を replay する。
        let mut transcript = Transcript::new(b"spartan-nizk-test");
        assert!(instance.verify(&io, &proof, &mut transcript));
        assert_eq!(
            transcript.challenge_field::<F>(b"probe"),
            prover_transcript.challenge_field::<F>(b"probe")
        );
    }

    #[test]
    fn proofs_are_randomized_but_all_verify() {
        let instance = square_mul_instance();
        let (io, witness) = example_io_witness();

        let prove = |seed: u64| {
            let mut transcript = Transcript::new(b"spartan-nizk-test");
            instance.prove(
                &io,
                &witness,
                &mut transcript,
                &mut StdRng::seed_from_u64(seed),
            )
        };
        let (lhs, rhs) = (prove(1), prove(2));

        // 乱数が違えば、witness を隠す commitment はどれも一致しない。
        assert_ne!(lhs.witness_com, rhs.witness_com);
        assert_ne!(lhs.va_com, rhs.va_com);
        assert_ne!(lhs.sc1.poly_coms, rhs.sc1.poly_coms);
        assert_ne!(lhs.sc2.eval_coms, rhs.sc2.eval_coms);
        assert_ne!(
            lhs.witness_opening.result_com,
            rhs.witness_opening.result_com
        );

        for proof in [lhs, rhs] {
            let mut transcript = Transcript::new(b"spartan-nizk-test");
            assert!(instance.verify(&io, &proof, &mut transcript));
        }
    }

    #[test]
    fn verifier_rejects_changed_statement_and_proof_values() {
        let instance = square_mul_instance();
        let (io, witness) = example_io_witness();

        let prove = || {
            let mut transcript = Transcript::new(b"spartan-nizk-test");
            instance.prove(
                &io,
                &witness,
                &mut transcript,
                &mut StdRng::seed_from_u64(42),
            )
        };
        let verify = |io: &[F], proof| {
            let mut transcript = Transcript::new(b"spartan-nizk-test");
            instance.verify(io, &proof, &mut transcript)
        };
        let g = instance.pcs.scalar.generator;

        assert!(!verify(&[F::from(4), F::from(75)], prove()));

        let mut proof = prove();
        proof.va_com += g;
        assert!(!verify(&io, proof));

        let mut proof = prove();
        proof.prod_com += g;
        assert!(!verify(&io, proof));

        let mut proof = prove();
        proof.sc1.eval_coms[0] += g;
        assert!(!verify(&io, proof));

        let mut proof = prove();
        proof.sc2.poly_coms[1] += g;
        assert!(!verify(&io, proof));

        let mut proof = prove();
        proof.witness_opening.result_com += g;
        assert!(!verify(&io, proof));

        let mut proof = prove();
        proof.sc2_equality.z += F::from(1);
        assert!(!verify(&io, proof));

        let mut proof = prove();
        proof.sc2.proofs.pop();
        assert!(!verify(&io, proof));
    }

    #[test]
    fn single_constraint_instance_has_a_zero_round_first_sumcheck() {
        // witness * 1 = out。制約 1 行なので constraint_vars = 0 で 1 本目の sumcheck は
        // ラウンドを持たず、最終 claim の commitment は初期値の単位元のまま。
        // half_len = 2, half_vars = 1 なので Hyrax の行列は 1 行 2 列になる。
        let r1cs = R1CS::<F, _>::new(
            DenseMatrix::from_usize([[0, 0, 1]]),
            DenseMatrix::from_usize([[1, 0, 0]]),
            DenseMatrix::from_usize([[0, 1, 0]]),
            R1CSStructure::new(1, 1, 1),
        );
        let instance =
            SpartanNizkInstance::<G1Projective, _>::encode::<G1Hasher>(SpartanR1CS::from(&r1cs))
                .unwrap();
        assert_eq!(instance.r1cs.constraint_vars(), 0);
        assert_eq!(instance.pcs.num_rows(), 1);
        assert_eq!(instance.pcs.rows.len(), 2);

        let io = [F::from(9)];
        let witness = [F::from(9)];
        let mut prover_transcript = Transcript::new(b"spartan-nizk-odd-test");
        let proof = instance.prove(
            &io,
            &witness,
            &mut prover_transcript,
            &mut StdRng::seed_from_u64(7),
        );
        assert!(proof.sc1.poly_coms.is_empty());

        let mut verifier_transcript = Transcript::new(b"spartan-nizk-odd-test");
        assert!(instance.verify(&io, &proof, &mut verifier_transcript));
    }

    #[test]
    #[should_panic(expected = "assignment does not satisfy the R1CS")]
    fn prove_rejects_unsatisfying_assignments() {
        let instance = square_mul_instance();
        let mut transcript = Transcript::new(b"spartan-nizk-test");
        instance.prove(
            &[F::from(3), F::from(74)],
            &[F::from(5), F::from(15)],
            &mut transcript,
            &mut StdRng::seed_from_u64(42),
        );
    }
}
