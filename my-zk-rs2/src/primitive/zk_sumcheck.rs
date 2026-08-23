//! Zero-knowledge sum-check protocol (Spartan 論文 8 章の構成)。
//!
//! 平文の sum-check ([`super::sumcheck`]) では round polynomial $s_j$ を
//! そのまま送るため、$g$ (ひいては witness) の情報が漏れる。zk 版では
//! $s_j$ の係数ベクトル $\vec{c}_j = (c_0, \ldots, c_d)$ を Pedersen commitment
//! で送り、verifier の 2 本の検査
//!
//! - $s_j(0) + s_j(1) = e_{j-1}$ (前ラウンドまでのクレーム)
//! - $s_j(r_j) = e_j$ (次ラウンドのクレーム)
//!
//! を係数ベクトルとの内積の形に直す:
//!
//! - $s_j(0) + s_j(1) = \langle \vec{c}_j, (2, 1, \ldots, 1) \rangle$
//! - $s_j(r_j) = \langle \vec{c}_j, (1, r_j, r_j^2, \ldots, r_j^d) \rangle$
//!
//! さらに verifier のランダムな重み $(w_0, w_1)$ で 2 本を 1 本に畳み、
//!
//! $w_0 e_{j-1} + w_1 e_j = \langle \vec{c}_j,\; w_0 (2, 1, \ldots, 1) + w_1 (1, r_j, \ldots, r_j^d) \rangle$
//!
//! をラウンドあたり 1 回の [`DotProductProof`] で証明する。左辺の commitment は
//! $w_0 C_{e_{j-1}} + w_1 C_{e_j}$ と準同型に計算できる。
//!
//! Hyrax は sum-check 全体を 1 回の dot product proof に乗せるが、Spartan は
//! 「ラウンドごとに独立した dot-product proof」に簡略化しており (論文 8 章)、
//! ここでもそれに従う。クレーム $e_j$ は commitment のまま伝搬し、最終クレーム
//! $g(r)$ の commitment を返す。その値の正しさの確認は呼び出し側の責務。

use ark_ec::CurveGroup;
use ark_ff::{Field, Zero};
use ark_std::{UniformRand, rand::Rng};

use super::pedersen::Pedersen;
use super::poly::MvPolynomial;
use super::sigma::DotProductProof;
use super::sumcheck::round_polynomial;
use super::transcript::Transcript;

/// zk sum-check の 1 ラウンド分のメッセージ。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ZkSumCheckRound<G: CurveGroup> {
    /// $s_j$ の係数ベクトルへの commitment $C_{s_j}$。
    pub comm_round_poly: G,
    /// 次ラウンドのクレーム $e_j = s_j(r_j)$ への commitment $C_{e_j}$。
    pub comm_eval: G,
    /// 2 本の検査を重み $(w_0, w_1)$ で畳んだ
    /// $w_0 e_{j-1} + w_1 e_j = \langle \vec{c}_j, w_0 \vec{a}_{sum} + w_1 \vec{a}_{eval} \rangle$ の証明。
    pub proof: DotProductProof<G>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ZkSumCheckProof<G: CurveGroup> {
    pub rounds: Vec<ZkSumCheckRound<G>>,
}

/// $s(0) + s(1) = 2 c_0 + c_1 + \cdots + c_d$ の重みベクトル $(2, 1, \ldots, 1)$。
fn sum_weights<F: Field>(len: usize) -> Vec<F> {
    (0..len)
        .map(|i| if i == 0 { F::from(2u64) } else { F::one() })
        .collect()
}

/// $s(r) = \langle \vec{c}, (1, r, r^2, \ldots) \rangle$ の重みベクトル。
fn eval_weights<F: Field>(r: &F, len: usize) -> Vec<F> {
    let mut powers = Vec::with_capacity(len);
    let mut power = F::one();
    for _ in 0..len {
        powers.push(power);
        power *= r;
    }
    powers
}

/// 2 本の検査を畳んだ重みベクトル $w_0 \vec{a}_{sum} + w_1 \vec{a}_{eval}$。
fn folded_weights<F: Field>(w: &[F; 2], r: &F, len: usize) -> Vec<F> {
    sum_weights::<F>(len)
        .iter()
        .zip(eval_weights(r, len))
        .map(|(sum, eval)| w[0] * sum + w[1] * eval)
        .collect()
}

/// 主張 $\sum_{x \in \{0,1\}^N} g(x) = e_0$ に対する zk sum-check prover。
///
/// $e_0$ は commitment $C_{e_0} = \mathrm{Com}(e_0; \rho_0)$ として扱い、
/// prover はその blind `claim_blind` $= \rho_0$ を知っている必要がある
/// (公開の主張なら $\rho_0 = 0$)。
///
/// 戻り値は (証明, 評価点 $r$, 最終クレーム $g(r)$, その blind)。
#[allow(clippy::type_complexity)]
pub fn prove_zk_sumcheck<G, const N: usize>(
    g: &MvPolynomial<G::ScalarField, N>,
    claim_blind: &G::ScalarField,
    degree_bound: usize,
    coeff_committer: &Pedersen<G>,
    scalar_committer: &Pedersen<G>,
    transcript: &mut Transcript,
    rng: &mut impl Rng,
) -> (
    ZkSumCheckProof<G>,
    [G::ScalarField; N],
    G::ScalarField,
    G::ScalarField,
)
where
    G: CurveGroup,
{
    transcript.append_usize(b"zk-sumcheck-num-rounds", N);
    transcript.append_usize(b"zk-sumcheck-degree-bound", degree_bound);

    let mut challenges = [G::ScalarField::zero(); N];
    let mut rounds = Vec::with_capacity(N);

    // (e_{j-1}, blind_{j-1})。最終ラウンド後は (g(r), その blind) になる。
    let mut eval = G::ScalarField::zero();
    let mut blind_eval = *claim_blind;

    for round in 0..N {
        let round_poly = round_polynomial(g, &challenges[..round], degree_bound);
        let coefficients = round_poly.coefficients(degree_bound + 1);

        let blind_poly = G::ScalarField::rand(rng);
        let comm_round_poly = coeff_committer.commit(&coefficients, &blind_poly);
        transcript.append_serializable(b"zk-sumcheck-round-poly-com", &comm_round_poly);

        let r: G::ScalarField = transcript.challenge_field(b"zk-sumcheck-challenge");
        challenges[round] = r;

        let next_eval = round_poly.eval(&[r]);
        let next_blind = G::ScalarField::rand(rng);
        let comm_eval = scalar_committer.commit_scalar(&next_eval, &next_blind);
        transcript.append_serializable(b"zk-sumcheck-eval-com", &comm_eval);

        // 2 本の検査
        //   s_j(0) + s_j(1) = e_{j-1},  s_j(r_j) = e_j
        // を重み (w_0, w_1) で 1 本の内積クレームに畳む。
        // 左辺 w_0 e_{j-1} + w_1 e_j の blind は w_0 ρ_{j-1} + w_1 ρ_j。
        let w: [G::ScalarField; 2] = transcript.challenge_field_array(b"zk-sumcheck-weight");
        let a = folded_weights(&w, &r, degree_bound + 1);
        let blind_target = w[0] * blind_eval + w[1] * next_blind;
        let proof = DotProductProof::prove(
            coeff_committer,
            scalar_committer,
            &coefficients,
            &blind_poly,
            &a,
            &blind_target,
            transcript,
            rng,
        );

        rounds.push(ZkSumCheckRound {
            comm_round_poly,
            comm_eval,
            proof,
        });
        eval = next_eval;
        blind_eval = next_blind;
    }

    (ZkSumCheckProof { rounds }, challenges, eval, blind_eval)
}

/// zk sum-check verifier。
///
/// クレーム $e_0$ の commitment `comm_claim` から出発し、各ラウンドの
/// dot-product proof を検証してクレームを縮小していく。
///
/// 成功時は評価点 $r$ と最終クレーム $g(r)$ の commitment を返す。
/// **その commitment が正しい値へのものかの確認は呼び出し側の責務**。
pub fn verify_zk_sumcheck<G, const N: usize>(
    comm_claim: &G,
    degree_bound: usize,
    proof: &ZkSumCheckProof<G>,
    coeff_committer: &Pedersen<G>,
    scalar_committer: &Pedersen<G>,
    transcript: &mut Transcript,
) -> Option<([G::ScalarField; N], G)>
where
    G: CurveGroup,
{
    if proof.rounds.len() != N {
        return None;
    }

    transcript.append_usize(b"zk-sumcheck-num-rounds", N);
    transcript.append_usize(b"zk-sumcheck-degree-bound", degree_bound);

    let mut challenges = [G::ScalarField::zero(); N];
    let mut comm_eval = *comm_claim;

    for (round, message) in proof.rounds.iter().enumerate() {
        transcript.append_serializable(b"zk-sumcheck-round-poly-com", &message.comm_round_poly);
        let r: G::ScalarField = transcript.challenge_field(b"zk-sumcheck-challenge");
        challenges[round] = r;
        transcript.append_serializable(b"zk-sumcheck-eval-com", &message.comm_eval);

        // w_0 e_{j-1} + w_1 e_j = <c_j, w_0 a_sum + w_1 a_eval>
        // 左辺の commitment は w_0 C_{e_{j-1}} + w_1 C_{e_j} と準同型に計算できる。
        let w: [G::ScalarField; 2] = transcript.challenge_field_array(b"zk-sumcheck-weight");
        let a = folded_weights(&w, &r, degree_bound + 1);
        let comm_target = comm_eval * w[0] + message.comm_eval * w[1];
        if !message.proof.verify(
            coeff_committer,
            scalar_committer,
            &a,
            &message.comm_round_poly,
            &comm_target,
            transcript,
        ) {
            return None;
        }

        comm_eval = message.comm_eval;
    }

    Some((challenges, comm_eval))
}

#[cfg(test)]
mod tests {
    use super::{prove_zk_sumcheck, verify_zk_sumcheck};
    use crate::primitive::{BoolHyperCube, Monomial, MvPolynomial, Pedersen, Transcript};
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use ark_std::test_rng;
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    fn f(x: u64) -> F {
        F::from(x)
    }

    fn committers() -> (Pedersen<G1Projective>, Pedersen<G1Projective>) {
        (
            Pedersen::setup::<G1Hasher>(b"zk-sumcheck-test", "coeff", 3).unwrap(),
            Pedersen::setup::<G1Hasher>(b"zk-sumcheck-test", "scalar", 1).unwrap(),
        )
    }

    fn example_polynomial() -> MvPolynomial<F, 3> {
        // g = 2*x0^2*x2 + 3*x1 + 7
        MvPolynomial::from_terms([
            (Monomial::new([2, 0, 1]), f(2)),
            (Monomial::new([0, 1, 0]), f(3)),
            (Monomial::new([0, 0, 0]), f(7)),
        ])
    }

    fn hypercube_sum(g: &MvPolynomial<F, 3>) -> F {
        BoolHyperCube::<3>::iter()
            .map(|point| g.eval(&point.to_field_point()))
            .sum()
    }

    #[test]
    fn zk_sumcheck_roundtrip_reduces_to_a_commitment_of_g_at_r() {
        let (coeff, scalar) = committers();
        let g = example_polynomial();
        let claim_blind = f(99);
        let comm_claim = scalar.commit_scalar(&hypercube_sum(&g), &claim_blind);
        let mut rng = test_rng();

        let mut prover_transcript = Transcript::new(b"zk-sumcheck-test");
        let (proof, r, final_eval, final_blind) = prove_zk_sumcheck(
            &g,
            &claim_blind,
            2,
            &coeff,
            &scalar,
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"zk-sumcheck-test");
        let (r_verifier, comm_final) = verify_zk_sumcheck::<_, 3>(
            &comm_claim,
            2,
            &proof,
            &coeff,
            &scalar,
            &mut verifier_transcript,
        )
        .expect("valid proof must verify");

        assert_eq!(r, r_verifier);
        assert_eq!(final_eval, g.eval(&r));
        assert_eq!(comm_final, scalar.commit_scalar(&final_eval, &final_blind));
    }

    #[test]
    fn zk_sumcheck_rejects_a_wrong_claim_commitment() {
        let (coeff, scalar) = committers();
        let g = example_polynomial();
        let claim_blind = f(99);
        let wrong_comm_claim = scalar.commit_scalar(&(hypercube_sum(&g) + f(1)), &claim_blind);
        let mut rng = test_rng();

        let mut prover_transcript = Transcript::new(b"zk-sumcheck-test");
        let (proof, _, _, _) = prove_zk_sumcheck(
            &g,
            &claim_blind,
            2,
            &coeff,
            &scalar,
            &mut prover_transcript,
            &mut rng,
        );

        let mut verifier_transcript = Transcript::new(b"zk-sumcheck-test");
        assert!(
            verify_zk_sumcheck::<_, 3>(
                &wrong_comm_claim,
                2,
                &proof,
                &coeff,
                &scalar,
                &mut verifier_transcript,
            )
            .is_none()
        );
    }

    #[test]
    fn zk_sumcheck_rejects_a_tampered_round() {
        let (coeff, scalar) = committers();
        let g = example_polynomial();
        let comm_claim = scalar.commit_scalar(&hypercube_sum(&g), &F::from(0));
        let mut rng = test_rng();

        let mut prover_transcript = Transcript::new(b"zk-sumcheck-test");
        let (mut proof, _, _, _) = prove_zk_sumcheck(
            &g,
            &F::from(0),
            2,
            &coeff,
            &scalar,
            &mut prover_transcript,
            &mut rng,
        );
        proof.rounds[1].comm_eval += scalar.generators[0];

        let mut verifier_transcript = Transcript::new(b"zk-sumcheck-test");
        assert!(
            verify_zk_sumcheck::<_, 3>(
                &comm_claim,
                2,
                &proof,
                &coeff,
                &scalar,
                &mut verifier_transcript,
            )
            .is_none()
        );
    }
}
