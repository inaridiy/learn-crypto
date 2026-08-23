use ark_ff::{Field, PrimeField};
use ark_serialize::CanonicalSerialize;

use super::poly::{MvPolynomial, lagrange};
use super::transcript::Transcript;

/// Sum-check protocol の証明。
///
/// `round_polys[j]` はラウンド $j$ の prover メッセージで、一変数多項式
///
/// $s_j(t) = \sum_{w \in \{0,1\}^{N-1-j}} g(r_0, \ldots, r_{j-1}, t, w)$
///
/// を表す ($r_0, \ldots, r_{j-1}$ はそれまでの verifier challenge)。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SumCheckProof<F: Field> {
    pub round_polys: Vec<MvPolynomial<F, 1>>,
}

/// 主張 $\sum_{x \in \{0,1\}^N} g(x) = \mathrm{claimed\_sum}$ に対する sum-check prover。
///
/// `degree_bound` は $g$ の各変数に関する次数の上界で、protocol parameter として
/// prover と verifier が共有する。各ラウンドで $s_j$ を transcript に書き込み、
/// challenge $r_j$ を得る (Fiat-Shamir で非対話化)。
///
/// 戻り値は証明と最終評価点 $r = (r_0, \ldots, r_{N-1})$。呼び出し側は $g(r)$ を
/// 別の手段 (commitment opening など) で検証者に提供する。
pub fn prove_sumcheck<F, const N: usize>(
    g: &MvPolynomial<F, N>,
    degree_bound: usize,
    transcript: &mut Transcript,
) -> (SumCheckProof<F>, [F; N])
where
    F: PrimeField + CanonicalSerialize,
{
    transcript.append_usize(b"sumcheck-num-rounds", N);
    transcript.append_usize(b"sumcheck-degree-bound", degree_bound);

    let mut challenges = [F::zero(); N];
    let mut round_polys = Vec::with_capacity(N);

    for round in 0..N {
        let round_poly = round_polynomial(g, &challenges[..round], degree_bound);
        transcript.append_polynomial(b"sumcheck-round-poly", &round_poly);
        challenges[round] = transcript.challenge_field(b"sumcheck-challenge");
        round_polys.push(round_poly);
    }

    (SumCheckProof { round_polys }, challenges)
}

/// ラウンド多項式 $s_j(t)$ を構成する。
///
/// $t = 0, 1, \ldots, d$ の各点で
/// $\sum_{w \in \{0,1\}^{N-j-1}} g(\mathrm{fixed}, t, w)$ を評価し、
/// $d + 1$ 点からの Lagrange 補間で一変数多項式に戻す。
pub(crate) fn round_polynomial<F: Field, const N: usize>(
    g: &MvPolynomial<F, N>,
    fixed: &[F],
    degree_bound: usize,
) -> MvPolynomial<F, 1> {
    let round = fixed.len();
    let num_suffix_vars = N - round - 1;

    let evaluations: Vec<(F, F)> = (0..=degree_bound)
        .map(|k| {
            let t = F::from(k as u64);
            let mut sum = F::zero();

            for suffix in 0..(1usize << num_suffix_vars) {
                let point: [F; N] = std::array::from_fn(|i| match i.cmp(&round) {
                    std::cmp::Ordering::Less => fixed[i],
                    std::cmp::Ordering::Equal => t,
                    std::cmp::Ordering::Greater => {
                        if (suffix >> (i - round - 1)) & 1 == 1 {
                            F::one()
                        } else {
                            F::zero()
                        }
                    }
                });
                sum += g.eval(&point);
            }

            (t, sum)
        })
        .collect();

    lagrange(&evaluations)
}

/// Sum-check verifier。
///
/// `claimed_sum` から出発し、各ラウンドで
///
/// - $\deg s_j \le d$
/// - $s_j(0) + s_j(1) = e_{j-1}$ (前ラウンドまでのクレーム)
///
/// を検査してクレームを $e_j = s_j(r_j)$ に縮小していく。
///
/// 成功時は評価点 $r$ と最終クレーム $g(r)$ を返す。**$g(r)$ がこの値になることの
/// 確認は呼び出し側の責務** (verifier は $g$ 自体を持たないため、oracle クエリや
/// commitment opening で確認する)。
pub fn verify_sumcheck<F, const N: usize>(
    claimed_sum: F,
    degree_bound: usize,
    proof: &SumCheckProof<F>,
    transcript: &mut Transcript,
) -> Option<([F; N], F)>
where
    F: PrimeField + CanonicalSerialize,
{
    if proof.round_polys.len() != N {
        return None;
    }

    transcript.append_usize(b"sumcheck-num-rounds", N);
    transcript.append_usize(b"sumcheck-degree-bound", degree_bound);

    let mut challenges = [F::zero(); N];
    let mut expected = claimed_sum;

    for (round, round_poly) in proof.round_polys.iter().enumerate() {
        if round_poly.degree().unwrap_or(0) > degree_bound {
            return None;
        }
        if round_poly.eval(&[F::zero()]) + round_poly.eval(&[F::one()]) != expected {
            return None;
        }

        transcript.append_polynomial(b"sumcheck-round-poly", round_poly);
        challenges[round] = transcript.challenge_field(b"sumcheck-challenge");
        expected = round_poly.eval(&[challenges[round]]);
    }

    Some((challenges, expected))
}

#[cfg(test)]
mod tests {
    use super::{prove_sumcheck, verify_sumcheck};
    use crate::primitive::{BoolHyperCube, Monomial, MvPolynomial, Transcript};
    use ark_bls12_381::Fr as F;

    fn f(x: u64) -> F {
        F::from(x)
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
    fn sumcheck_roundtrip_accepts_and_reduces_to_g_at_r() {
        let g = example_polynomial();
        let claimed_sum = hypercube_sum(&g);

        let mut prover_transcript = Transcript::new(b"sumcheck-test");
        let (proof, r) = prove_sumcheck(&g, 2, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"sumcheck-test");
        let (r_verifier, final_claim) =
            verify_sumcheck::<F, 3>(claimed_sum, 2, &proof, &mut verifier_transcript)
                .expect("valid proof must verify");

        assert_eq!(r, r_verifier);
        assert_eq!(final_claim, g.eval(&r));
    }

    #[test]
    fn sumcheck_rejects_wrong_claimed_sum() {
        let g = example_polynomial();
        let claimed_sum = hypercube_sum(&g) + f(1);

        let mut prover_transcript = Transcript::new(b"sumcheck-test");
        let (proof, _) = prove_sumcheck(&g, 2, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"sumcheck-test");
        assert!(
            verify_sumcheck::<F, 3>(claimed_sum, 2, &proof, &mut verifier_transcript).is_none()
        );
    }

    #[test]
    fn sumcheck_rejects_tampered_round_polynomial() {
        let g = example_polynomial();
        let claimed_sum = hypercube_sum(&g);

        let mut prover_transcript = Transcript::new(b"sumcheck-test");
        let (mut proof, _) = prove_sumcheck(&g, 2, &mut prover_transcript);
        proof.round_polys[1] += MvPolynomial::constant(f(1));

        let mut verifier_transcript = Transcript::new(b"sumcheck-test");
        assert!(
            verify_sumcheck::<F, 3>(claimed_sum, 2, &proof, &mut verifier_transcript).is_none()
        );
    }

    #[test]
    fn sumcheck_rejects_round_polynomial_over_degree_bound() {
        let g = example_polynomial();
        let claimed_sum = hypercube_sum(&g);

        // degree bound 1 に対して g は x0 について次数 2 なので拒否される。
        let mut prover_transcript = Transcript::new(b"sumcheck-test");
        let (proof, _) = prove_sumcheck(&g, 2, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"sumcheck-test");
        assert!(
            verify_sumcheck::<F, 3>(claimed_sum, 1, &proof, &mut verifier_transcript).is_none()
        );
    }
}
