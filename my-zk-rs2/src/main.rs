#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

use ark_ec::{CurveGroup, hashing::HashToCurve, hashing::HashToCurveError};
use ark_ff::{Field, One, Zero};
use ark_serialize::CanonicalSerialize;
use ark_std::{UniformRand, rand::Rng};
use my_zk_rs::primitive::{
    BoolHyperCube, EqualityProof, KnowledgeProof, Matrix, MvPolynomial, Pedersen, ProductProof,
    R1CS, Transcript, ZkSumCheckProof,
    hyrax::{HyraxPCS, HyraxPCSCommitment, HyraxPCSProof},
    mle_from_hypercube_evaluations, mle_from_matrix, prove_zk_sumcheck, teq, verify_zk_sumcheck,
};

/// sum-check #1 の対象 $G_\tau$ の各変数についての次数上界。
/// $\bar{A} \cdot \bar{B}$ (multilinear の積で次数 2) と $\widetilde{eq}$ (次数 1) の積。
const FIRST_SUMCHECK_DEGREE: usize = 3;

/// sum-check #2 の対象 $M$ の各変数についての次数上界。multilinear 2 つの積。
const SECOND_SUMCHECK_DEGREE: usize = 2;

#[derive(Clone, Debug)]
pub struct SpartanNIZKGens<G: CurveGroup, const HW: usize>
where
    [(); 1 << HW]:,
{
    pub pcs: HyraxPCS<G, HW>,
    /// blind 生成元は `pcs` と共有される。
    pub round_poly_committer: Pedersen<G>,
}

impl<G: CurveGroup, const HW: usize> SpartanNIZKGens<G, HW>
where
    [(); 1 << HW]:,
    [(); 1 << (HW * 2)]:,
{
    pub fn setup<H: HashToCurve<G>>(domain: &[u8]) -> Result<Self, HashToCurveError> {
        Ok(Self {
            pcs: HyraxPCS::setup::<H>(domain)?,
            round_poly_committer: Pedersen::setup::<H>(
                domain,
                "spartan-round-poly",
                FIRST_SUMCHECK_DEGREE + 1,
            )?,
        })
    }
}

/// SpartanNIZK の証明。witness に依存する値はすべて commitment。
#[derive(Clone, Debug)]
pub struct SpartanNIZKProof<G: CurveGroup, const HW: usize>
where
    [(); 1 << HW]:,
{
    /// $\tilde{w}$ への Hyrax commitment。
    pub witness_com: HyraxPCSCommitment<G, HW>,
    /// zk sum-check #1: $\sum_x G_\tau(x) = 0$。
    pub zk_sumcheck1: ZkSumCheckProof<G>,
    /// $C_{v_A}, C_{v_B}, C_{v_C}$。
    pub comm_v_a: G,
    pub comm_v_b: G,
    pub comm_v_c: G,
    /// $C_{v_{AB}}$ ($v_{AB} = v_A \cdot v_B$)。
    pub comm_v_ab: G,
    /// $C_{v_C}$ の開示値の proof of knowledge。
    pub pok_v_c: KnowledgeProof<G>,
    /// $v_{AB} = v_A \cdot v_B$ の product proof。
    pub product_proof: ProductProof<G>,
    /// $e_x = (v_{AB} - v_C) \cdot \widetilde{eq}(r_x, \tau)$ の equality proof。
    pub eq1_proof: EqualityProof<G>,
    /// zk sum-check #2: $\sum_y M(y) = r_A v_A + r_B v_B + r_C v_C$。
    pub zk_sumcheck2: ZkSumCheckProof<G>,
    /// $\tilde{w}(r_y')$ の Hyrax opening proof。評価値は `result_com` に隠れたまま。
    pub witness_opening: HyraxPCSProof<G>,
    /// $e_y = m \cdot \tilde{Z}(r_y)$ の equality proof。
    pub eq2_proof: EqualityProof<G>,
}

/// $\bar{A}(x) = \sum_{y \in \{0,1\}^M} \tilde{A}(x, y) \cdot \tilde{Z}(y)$
///
/// ブール点上では $\tilde{Z}(y) = z_y$ なので、$\bar{A}$ は $Az$ の
/// multilinear extension に一致する。
fn calc_bar_matrix<F, const N: usize, const M: usize>(
    matrix: &Matrix<F, N, M>,
    tassignment: &MvPolynomial<F, M>,
) -> MvPolynomial<F, N>
where
    F: Field,
    [(); 1 << N]:,
    [(); 1 << M]:,
    [(); 1 << (N + M)]:,
{
    let matrix_mle = mle_from_matrix(matrix);
    let mut result = MvPolynomial::zero();

    for y in BoolHyperCube::<M>::iter() {
        let y = y.to_field_point();
        let assignment_value = tassignment.eval(&y);
        if assignment_value.is_zero() {
            continue;
        }

        result += matrix_mle.curry_suffix(&y).scale(assignment_value);
    }

    result
}

/// 行列 MLE の row 変数を固定した一部評価 $\tilde{A}(r_x, \cdot)$ を返す。
///
/// Prover は sum-check #2 の対象多項式の構成に、verifier は
/// $\tilde{A}(r_x, r_y)$ の評価 (非簡潔な $O(n)$ の仕事) に使う。
fn matrix_row_mle<F, const N: usize, const M: usize>(
    matrix: &Matrix<F, N, M>,
    row_point: &[F; N],
) -> MvPolynomial<F, M>
where
    F: Field,
    [(); 1 << N]:,
    [(); 1 << M]:,
    [(); 1 << (N + M)]:,
{
    mle_from_matrix(matrix).curry_prefix(row_point)
}

/// 証明対象の statement (R1CS instance と io) を transcript に束縛する。
fn init_transcript<F, const S: usize, const M: usize>(r1cs: &R1CS<F, S, M>, io: &[F]) -> Transcript
where
    F: Field + CanonicalSerialize,
    [(); 1 << S]:,
    [(); 1 << M]:,
{
    let mut transcript = Transcript::new(b"spartan/nizk");
    transcript.append_bytes(b"field", std::any::type_name::<F>().as_bytes());
    transcript.append_usize(b"constraint-bits", S);
    transcript.append_usize(b"var-bits", M);
    transcript.append_usize(b"num-constraints", r1cs.structure.num_constraints);
    transcript.append_usize(b"num-io", r1cs.structure.num_io);
    transcript.append_usize(b"num-witness", r1cs.structure.num_witness);
    transcript.append_matrix(b"r1cs-a", &r1cs.a);
    transcript.append_matrix(b"r1cs-b", &r1cs.b);
    transcript.append_matrix(b"r1cs-c", &r1cs.c);
    for value in io {
        transcript.append_field(b"io", value);
    }
    transcript
}

pub fn prove<G, const S: usize, const HW: usize>(
    gens: &SpartanNIZKGens<G, HW>,
    r1cs: &R1CS<G::ScalarField, S, { HW * 2 + 1 }>,
    io: &[G::ScalarField],
    witness: &[G::ScalarField],
    rng: &mut impl Rng,
) -> SpartanNIZKProof<G, HW>
where
    G: CurveGroup,
    [(); 1 << S]:,
    [(); 1 << HW]:,
    [(); HW * 2]:,
    [(); 1 << (HW * 2)]:,
    [(); 1 << (HW * 2 + 1)]:,
    [(); 1 << (S + (HW * 2 + 1))]:,
{
    let z = r1cs.assignment(io, witness);
    assert!(r1cs.is_sat(&z), "the witness does not satisfy the R1CS");

    let mut transcript = init_transcript(r1cs, io);
    let scalar_committer = &gens.pcs.scalar;

    // 1. witness の MLE w̃ (= z の後半の MLE) に Hyrax でコミットする。
    let half = 1 << (HW * 2);
    let witness_table: [G::ScalarField; 1 << (HW * 2)] = std::array::from_fn(|i| z[half + i]);
    let witness_blinds: [G::ScalarField; 1 << HW] =
        std::array::from_fn(|_| G::ScalarField::rand(rng));
    let witness_com = gens.pcs.commit(&witness_table, witness_blinds);
    for commitment in &witness_com {
        transcript.append_serializable(b"witness-commitment", commitment);
    }

    // 2. $\tau \in \mathbb{F}^S$
    let tau: [G::ScalarField; S] = transcript.challenge_field_array(b"tau");

    // 3. zk sum-check #1: $\sum_x (\bar{A} \bar{B} - \bar{C})(x) \cdot \widetilde{eq}(x, \tau) = 0$
    //    クレーム 0 は公開値なので commitment は $\mathrm{Com}(0; 0)$ (= 単位元)。
    let z_mle = mle_from_hypercube_evaluations::<G::ScalarField, { HW * 2 + 1 }>(z);
    let bar_a = calc_bar_matrix(&r1cs.a, &z_mle);
    let bar_b = calc_bar_matrix(&r1cs.b, &z_mle);
    let bar_c = calc_bar_matrix(&r1cs.c, &z_mle);
    let g1 = (&bar_a * &bar_b - &bar_c) * teq(&tau);
    let (zk_sumcheck1, r_x, e_x, e_x_blind) = prove_zk_sumcheck(
        &g1,
        &G::ScalarField::zero(),
        FIRST_SUMCHECK_DEGREE,
        &gens.round_poly_committer,
        scalar_committer,
        &mut transcript,
        rng,
    );

    // 4. $v_A = \bar{A}(r_x), v_B = \bar{B}(r_x), v_C = \bar{C}(r_x)$ への commitment。
    let v_a = bar_a.eval(&r_x);
    let v_b = bar_b.eval(&r_x);
    let v_c = bar_c.eval(&r_x);
    let blind_v_a = G::ScalarField::rand(rng);
    let blind_v_b = G::ScalarField::rand(rng);
    let blind_v_c = G::ScalarField::rand(rng);
    let blind_v_ab = G::ScalarField::rand(rng);
    let comm_v_a = scalar_committer.commit_scalar(&v_a, &blind_v_a);
    let comm_v_b = scalar_committer.commit_scalar(&v_b, &blind_v_b);
    let comm_v_c = scalar_committer.commit_scalar(&v_c, &blind_v_c);
    let comm_v_ab = scalar_committer.commit_scalar(&(v_a * v_b), &blind_v_ab);
    transcript.append_serializable(b"comm-v-a", &comm_v_a);
    transcript.append_serializable(b"comm-v-b", &comm_v_b);
    transcript.append_serializable(b"comm-v-c", &comm_v_c);

    // $C_{v_C}$ の開示値の proof of knowledge
    // $v_A, v_B, v_{AB}$ は product proof が knowledge の証明を兼ねる。
    let pok_v_c = KnowledgeProof::prove(scalar_committer, &v_c, &blind_v_c, &mut transcript, rng);

    // $v_{AB} = v_A \cdot v_B$ の product proof。
    let product_proof = ProductProof::prove(
        scalar_committer,
        &v_a,
        &blind_v_a,
        &v_b,
        &blind_v_b,
        &blind_v_ab,
        &mut transcript,
        rng,
    );

    // $e_x = (v_{AB} - v_C) \cdot \widetilde{eq}(r_x, \tau)$ の equality proof。
    // 右辺の commitment は $(C_{v_{AB}} - C_{v_C}) \cdot \widetilde{eq}(r_x, \tau)$ と
    // 準同型に作れる。
    let eq_tau_r_x = teq(&tau).eval(&r_x);
    let comm_e_x = scalar_committer.commit_scalar(&e_x, &e_x_blind);
    let comm_rhs1 = (comm_v_ab - comm_v_c) * eq_tau_r_x;
    let blind_rhs1 = (blind_v_ab - blind_v_c) * eq_tau_r_x;
    let eq1_proof = EqualityProof::prove(
        &gens.pcs.scalar.blind,
        &comm_e_x,
        &comm_rhs1,
        &(e_x_blind - blind_rhs1),
        &mut transcript,
        rng,
    );

    // 5. $r_A, r_B, r_C$
    let r_a: G::ScalarField = transcript.challenge_field(b"r_a");
    let r_b: G::ScalarField = transcript.challenge_field(b"r_b");
    let r_c: G::ScalarField = transcript.challenge_field(b"r_c");

    // 6. zk sum-check #2:
    //    $\sum_y (r_A \tilde{A}(r_x, y) + r_B \tilde{B}(r_x, y) + r_C \tilde{C}(r_x, y)) \cdot \tilde{Z}(y)
    //      = r_A v_A + r_B v_B + r_C v_C$
    //    クレームの blind は $r_A \rho_{v_A} + r_B \rho_{v_B} + r_C \rho_{v_C}$。
    let a_r_x = matrix_row_mle(&r1cs.a, &r_x);
    let b_r_x = matrix_row_mle(&r1cs.b, &r_x);
    let c_r_x = matrix_row_mle(&r1cs.c, &r_x);
    let g2 =
        (a_r_x.clone().scale(r_a) + b_r_x.clone().scale(r_b) + c_r_x.clone().scale(r_c)) * &z_mle;
    let claim2_blind = r_a * blind_v_a + r_b * blind_v_b + r_c * blind_v_c;
    let (zk_sumcheck2, r_y, e_y, e_y_blind) = prove_zk_sumcheck(
        &g2,
        &claim2_blind,
        SECOND_SUMCHECK_DEGREE,
        &gens.round_poly_committer,
        scalar_committer,
        &mut transcript,
        rng,
    );

    // 7. $\tilde{w}(r_y')$ の Hyrax opening proof ($r_y'$ = $r_y$ の下位 HW*2 変数)。
    //    評価値は `result_com` に blind 付きで隠れたまま。
    let r_y_witness: [G::ScalarField; HW * 2] = std::array::from_fn(|i| r_y[i]);
    let result_blind = G::ScalarField::rand(rng);
    let witness_opening = gens.pcs.prove_with_transcript(
        &witness_com,
        &witness_table,
        &witness_blinds,
        &r_y_witness,
        &result_blind,
        &mut transcript,
        rng,
    );

    // 8. $e_y = m \cdot \tilde{Z}(r_y)$ の equality proof。
    //    $m = r_A \tilde{A}(r_x, r_y) + r_B \tilde{B}(r_x, r_y) + r_C \tilde{C}(r_x, r_y)$ と
    //    $\tilde{Z}(r_y)$ の公開部分は誰でも計算でき、witness 部分は `result_com`。
    let m = r_a * a_r_x.eval(&r_y) + r_b * b_r_x.eval(&r_y) + r_c * c_r_x.eval(&r_y);
    let r_y_last = r_y[HW * 2];
    let mut public_table = [G::ScalarField::zero(); 1 << (HW * 2)];
    public_table[..io.len()].copy_from_slice(io);
    public_table[io.len()] = G::ScalarField::one();
    let public_eval = mle_from_hypercube_evaluations::<G::ScalarField, { HW * 2 }>(public_table)
        .eval(&r_y_witness);
    // $\mathrm{Com}(m \tilde{Z}(r_y)) = m \cdot (\mathrm{Com}((1 - r_{y,last}) \widetilde{(io,1)}(r_y'); 0) + r_{y,last} \cdot \mathrm{result\_com})$
    let comm_z_eval = scalar_committer.commit_scalar(
        &((G::ScalarField::one() - r_y_last) * public_eval),
        &G::ScalarField::zero(),
    ) + witness_opening.result_com * r_y_last;
    let comm_rhs2 = comm_z_eval * m;
    let blind_rhs2 = m * r_y_last * result_blind;
    let comm_e_y = scalar_committer.commit_scalar(&e_y, &e_y_blind);
    let eq2_proof = EqualityProof::prove(
        &gens.pcs.scalar.blind,
        &comm_e_y,
        &comm_rhs2,
        &(e_y_blind - blind_rhs2),
        &mut transcript,
        rng,
    );

    SpartanNIZKProof {
        witness_com,
        zk_sumcheck1,
        comm_v_a,
        comm_v_b,
        comm_v_c,
        comm_v_ab,
        pok_v_c,
        product_proof,
        eq1_proof,
        zk_sumcheck2,
        witness_opening,
        eq2_proof,
    }
}

pub fn verify<G, const S: usize, const HW: usize>(
    gens: &SpartanNIZKGens<G, HW>,
    r1cs: &R1CS<G::ScalarField, S, { HW * 2 + 1 }>,
    io: &[G::ScalarField],
    proof: &SpartanNIZKProof<G, HW>,
) -> bool
where
    G: CurveGroup,
    [(); 1 << S]:,
    [(); 1 << HW]:,
    [(); HW * 2]:,
    [(); 1 << (HW * 2)]:,
    [(); 1 << (HW * 2 + 1)]:,
    [(); 1 << (S + (HW * 2 + 1))]:,
{
    if io.len() != r1cs.structure.num_io {
        return false;
    }

    let mut transcript = init_transcript(r1cs, io);
    let scalar_committer = &gens.pcs.scalar;

    for commitment in &proof.witness_com {
        transcript.append_serializable(b"witness-commitment", commitment);
    }
    let tau: [G::ScalarField; S] = transcript.challenge_field_array(b"tau");

    // zk sum-check #1。クレーム 0 の commitment は $\mathrm{Com}(0; 0)$ = 単位元。
    let comm_claim1 = G::zero();
    let Some((r_x, comm_e_x)) = verify_zk_sumcheck::<_, S>(
        &comm_claim1,
        FIRST_SUMCHECK_DEGREE,
        &proof.zk_sumcheck1,
        &gens.round_poly_committer,
        scalar_committer,
        &mut transcript,
    ) else {
        return false;
    };

    transcript.append_serializable(b"comm-v-a", &proof.comm_v_a);
    transcript.append_serializable(b"comm-v-b", &proof.comm_v_b);
    transcript.append_serializable(b"comm-v-c", &proof.comm_v_c);

    // $C_{v_C}$ の開示値の proof of knowledge。
    if !proof
        .pok_v_c
        .verify(scalar_committer, &proof.comm_v_c, &mut transcript)
    {
        return false;
    }

    // $v_{AB} = v_A \cdot v_B$
    if !proof.product_proof.verify(
        scalar_committer,
        &proof.comm_v_a,
        &proof.comm_v_b,
        &proof.comm_v_ab,
        &mut transcript,
    ) {
        return false;
    }

    // $e_x = (v_{AB} - v_C) \cdot \widetilde{eq}(r_x, \tau)$
    let eq_tau_r_x = teq(&tau).eval(&r_x);
    let comm_rhs1 = (proof.comm_v_ab - proof.comm_v_c) * eq_tau_r_x;
    if !proof.eq1_proof.verify(
        &gens.pcs.scalar.blind,
        &comm_e_x,
        &comm_rhs1,
        &mut transcript,
    ) {
        return false;
    }

    let r_a: G::ScalarField = transcript.challenge_field(b"r_a");
    let r_b: G::ScalarField = transcript.challenge_field(b"r_b");
    let r_c: G::ScalarField = transcript.challenge_field(b"r_c");

    // zk sum-check #2。クレームの commitment は準同型に計算できる。
    let comm_claim2 = proof.comm_v_a * r_a + proof.comm_v_b * r_b + proof.comm_v_c * r_c;
    let Some((r_y, comm_e_y)) = verify_zk_sumcheck::<_, { HW * 2 + 1 }>(
        &comm_claim2,
        SECOND_SUMCHECK_DEGREE,
        &proof.zk_sumcheck2,
        &gens.round_poly_committer,
        scalar_committer,
        &mut transcript,
    ) else {
        return false;
    };

    // $\tilde{w}(r_y')$ の opening proof を検証する。評価値は `result_com` の中。
    let r_y_witness: [G::ScalarField; HW * 2] = std::array::from_fn(|i| r_y[i]);
    if !gens.pcs.verify_with_transcript(
        &proof.witness_com,
        &r_y_witness,
        &proof.witness_opening,
        &mut transcript,
    ) {
        return false;
    }

    // 非簡潔な部分: verifier が $\tilde{A}(r_x, r_y)$ などを行列の MLE から自力で評価する。
    let m = r_a * matrix_row_mle(&r1cs.a, &r_x).eval(&r_y)
        + r_b * matrix_row_mle(&r1cs.b, &r_x).eval(&r_y)
        + r_c * matrix_row_mle(&r1cs.c, &r_x).eval(&r_y);

    // $\tilde{Z}(r_y) = (1 - r_{y,last}) \cdot \widetilde{(io,1)}(r_y') + r_{y,last} \cdot \tilde{w}(r_y')$
    let r_y_last = r_y[HW * 2];
    let mut public_table = [G::ScalarField::zero(); 1 << (HW * 2)];
    public_table[..io.len()].copy_from_slice(io);
    public_table[io.len()] = G::ScalarField::one();
    let public_eval = mle_from_hypercube_evaluations::<G::ScalarField, { HW * 2 }>(public_table)
        .eval(&r_y_witness);
    let comm_z_eval = scalar_committer.commit_scalar(
        &((G::ScalarField::one() - r_y_last) * public_eval),
        &G::ScalarField::zero(),
    ) + proof.witness_opening.result_com * r_y_last;
    let comm_rhs2 = comm_z_eval * m;

    // $e_y = m \cdot \tilde{Z}(r_y)$
    proof.eq2_proof.verify(
        &gens.pcs.scalar.blind,
        &comm_e_y,
        &comm_rhs2,
        &mut transcript,
    )
}

mod example {
    use super::*;
    use ark_bls12_381::Fr as F;
    use my_zk_rs::primitive::{Matrix, R1CSStructure};

    /// $out = x^3 + x + 5$ を表す R1CS。io = $(x, out)$、witness = $(i_1, i_2, i_3)$。
    ///
    /// $z = (x, out, 1, 0 \mid i_1, i_2, i_3, 0)$ に対する制約:
    ///
    /// - $x \cdot x = i_1$
    /// - $i_1 \cdot x = i_2$
    /// - $(x + i_2) \cdot 1 = i_3$
    /// - $(i_3 + 5) \cdot 1 = out$
    pub fn cubic_r1cs() -> R1CS<F, 2, 3> {
        R1CS::new(
            Matrix::from_usize([
                [1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 5, 0, 0, 0, 1, 0],
            ]),
            Matrix::from_usize([
                [1, 0, 0, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0],
            ]),
            Matrix::from_usize([
                [0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
            ]),
            R1CSStructure {
                num_constraints: 4,
                num_io: 2,
                num_witness: 3,
            },
        )
    }
}

pub fn main() {
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    let gens = SpartanNIZKGens::<G1Projective, 1>::setup::<G1Hasher>(b"spartan-nizk").unwrap();
    let r1cs = example::cubic_r1cs();

    // x = 3 に対して out = 3^3 + 3 + 5 = 35。
    let io = [3, 35].map(F::from);
    let witness = [9, 27, 30].map(F::from);

    let mut rng = ark_std::test_rng();
    let proof = prove(&gens, &r1cs, &io, &witness, &mut rng);
    let verified = verify(&gens, &r1cs, &io, &proof);

    println!("SpartanNIZK verified: {verified}");
    assert!(verified);
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use ark_std::test_rng;
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    fn example_setup() -> (
        SpartanNIZKGens<G1Projective, 1>,
        R1CS<F, 2, 3>,
        [F; 2],
        [F; 3],
    ) {
        let gens = SpartanNIZKGens::<G1Projective, 1>::setup::<G1Hasher>(b"spartan-nizk").unwrap();
        let r1cs = example::cubic_r1cs();
        let io = [3, 35].map(F::from);
        let witness = [9, 27, 30].map(F::from);
        (gens, r1cs, io, witness)
    }

    #[test]
    fn spartan_nizk_accepts_valid_proof() {
        let (gens, r1cs, io, witness) = example_setup();
        let mut rng = test_rng();

        let proof = prove(&gens, &r1cs, &io, &witness, &mut rng);

        assert!(verify(&gens, &r1cs, &io, &proof));
    }

    #[test]
    fn spartan_nizk_rejects_different_io() {
        let (gens, r1cs, io, witness) = example_setup();
        let mut rng = test_rng();

        let proof = prove(&gens, &r1cs, &io, &witness, &mut rng);
        let other_io = [3, 36].map(F::from);

        assert!(!verify(&gens, &r1cs, &other_io, &proof));
    }

    #[test]
    fn spartan_nizk_rejects_tampered_value_commitment() {
        let (gens, r1cs, io, witness) = example_setup();
        let mut rng = test_rng();

        let mut proof = prove(&gens, &r1cs, &io, &witness, &mut rng);
        proof.comm_v_a += gens.pcs.scalar.generators[0];

        assert!(!verify(&gens, &r1cs, &io, &proof));
    }

    #[test]
    fn spartan_nizk_rejects_tampered_witness_opening() {
        let (gens, r1cs, io, witness) = example_setup();
        let mut rng = test_rng();

        let mut proof = prove(&gens, &r1cs, &io, &witness, &mut rng);
        proof.witness_opening.result_com += gens.pcs.scalar.generators[0];

        assert!(!verify(&gens, &r1cs, &io, &proof));
    }

    #[test]
    fn spartan_nizk_rejects_tampered_equality_proof() {
        let (gens, r1cs, io, witness) = example_setup();
        let mut rng = test_rng();

        let mut proof = prove(&gens, &r1cs, &io, &witness, &mut rng);
        proof.eq2_proof.z += F::from(1);

        assert!(!verify(&gens, &r1cs, &io, &proof));
    }

    #[test]
    fn spartan_nizk_rejects_tampered_knowledge_proof() {
        let (gens, r1cs, io, witness) = example_setup();
        let mut rng = test_rng();

        let mut proof = prove(&gens, &r1cs, &io, &witness, &mut rng);
        proof.pok_v_c.z1 += F::from(1);

        assert!(!verify(&gens, &r1cs, &io, &proof));
    }

    #[test]
    #[should_panic(expected = "does not satisfy")]
    fn spartan_nizk_prover_rejects_bad_witness() {
        let (gens, r1cs, io, _) = example_setup();
        let mut rng = test_rng();

        let bad_witness = [9, 27, 31].map(F::from);
        prove(&gens, &r1cs, &io, &bad_witness, &mut rng);
    }
}
