use std::collections::BTreeSet;

use ark_ff::PrimeField;

use crate::primitive::{EvalsUniPoly, MultilinearPoly, Transcript, UniPoly, hypercube_size};

pub struct SumcheckProof<F: PrimeField> {
    polys: Vec<EvalsUniPoly<F>>,
}

/// `prove` の出力。proof に加えて、prover 側だけが持つ終点の情報を返す。
pub struct SumcheckOutput<F: PrimeField> {
    pub proof: SumcheckProof<F>,
    /// 各ラウンドの challenge `r = (r_0, ..., r_{n-1})`。sumcheck の終点。
    pub point: Vec<F>,
    /// 終点での各因子の値 `[p_0(r), ..., p_{k-1}(r)]`。
    pub final_evals: Vec<F>,
}

impl<F: PrimeField> SumcheckProof<F> {
    /// `polys` が sum-check の入力として整合しているか確かめ、ラウンド数(変数の数)を返す。
    pub fn assert_shape<P: MultilinearPoly<F>>(polys: &[P], degree: usize) -> usize {
        assert!(!polys.is_empty(), "sum-check needs at least one polynomial");
        assert!(degree >= 1, "round polynomial degree must be at least 1");

        let num_rounds = polys[0].vars();
        for (i, p) in polys.iter().enumerate().skip(1) {
            assert_eq!(p.vars(), num_rounds, "polynomial {i} dimension differs");
        }
        num_rounds
    }

    /// 現在の `polys` に対するラウンド多項式
    /// $s(t) = \sum_{x'} \mathrm{comb}(p_0(t, x'), \ldots, p_{k-1}(t, x'))$ を計算する。
    ///
    /// ラウンド多項式は `t = 0, 2, 3, ..., degree` の `degree` 点で評価し、
    /// `t = 1` は `claim - s(0)` から導く。`degree` が実際の次数より小さいと
    /// 補間が本物の多項式と一致せず、検証に失敗する。
    ///
    /// hypercube 上のペアのうち、どれかの因子が非零になるものだけを走る。
    /// 全因子が zero のペアは `comb(0, ..., 0)` を返すだけなので、個数を数えてまとめて足す。
    /// 1 ラウンドのコストは $O(k \cdot \mathrm{nnz} \cdot \text{degree})$。
    pub fn round_polynomial<P: MultilinearPoly<F>, Comb: Fn(&[F]) -> F>(
        polys: &[P],
        claim: F,
        degree: usize,
        comb: &Comb,
    ) -> EvalsUniPoly<F> {
        let zeros = vec![F::zero(); polys.len()];

        // evals[t] = s(t)
        let mut evals = vec![F::zero(); degree + 1];

        // どれかの因子が lo か hi で非零になるペア i = index >> 1 の集合。
        let active_pairs: BTreeSet<usize> = polys
            .iter()
            .flat_map(|p| p.nonzero_entries().map(|(index, _)| index >> 1))
            .collect();

        for &i in &active_pairs {
            let (lo, hi) = (2 * i, 2 * i + 1);
            let p_lo: Vec<F> = polys.iter().map(|p| p.eval_index(lo)).collect();
            let p_hi: Vec<F> = polys.iter().map(|p| p.eval_index(hi)).collect();
            let delta: Vec<F> = p_hi.iter().zip(&p_lo).map(|(hi, lo)| *hi - lo).collect();

            evals[0] += comb(&p_lo);
            // evals[1] は claim から導くのでここでは計算しない。

            // multilinear なので p(t) = p(0) + t * (p(1) - p(0))。
            for t in 2..=degree {
                let p_t: Vec<F> = p_lo
                    .iter()
                    .zip(&delta)
                    .map(|(lo, d)| *lo + F::from(t as u64) * d)
                    .collect();
                evals[t] += comb(&p_t);
            }
        }

        // 全因子が zero のペアは、どの t でも comb(0, ..., 0) を足すだけ。
        let zero_pairs = hypercube_size(polys[0].vars()) / 2 - active_pairs.len();
        let zero_term = F::from(zero_pairs as u64) * comb(&zeros);
        evals[0] += zero_term;
        for t in 2..=degree {
            evals[t] += zero_term;
        }

        // s(0) + s(1) = claim
        evals[1] = claim - evals[0];

        EvalsUniPoly::new(evals)
    }

    /// 不正な proof を panic せず拒否する verifier 向けの検証 API。
    ///
    /// `degree` はプロトコルが定めるラウンド多項式の次数上限で、verifier が知っている値。
    /// 成功時は `(最終 claim, r)` を返す。
    pub fn verify(
        &self,
        sum: F,
        degree: usize,
        transcript: &mut Transcript,
    ) -> Option<(F, Vec<F>)> {
        let mut h = sum;
        let mut r: Vec<F> = Vec::with_capacity(self.polys.len());

        for p in &self.polys {
            if p.degree() > degree {
                return None;
            }

            let p = &p.to_coeffs();
            let eval_zero = p.eval(F::zero());
            let eval_one = p.eval(F::one());
            if eval_one + eval_zero != h {
                return None;
            }

            transcript.append_serializable(b"round_poly", p);
            let r_i = transcript.challenge_field::<F>(b"challenge_r");
            r.push(r_i);

            h = p.eval(r_i);
        }

        Some((h, r))
    }

    /// `comb(p_0(x), ..., p_{k-1}(x))` を hypercube 上で足し上げる sum-check。
    ///
    /// `comb` は multilinear 多項式 `polys` の評価値を変数とする任意の多項式で、
    /// `degree` はその各変数についての次数(= ラウンド多項式の次数)。
    /// ラウンド多項式の計算は [`Self::round_polynomial`] を参照。
    pub fn prove<P: MultilinearPoly<F>, Comb: Fn(&[F]) -> F>(
        mut polys: Vec<P>,
        mut claim: F,
        degree: usize,
        comb: Comb,
        transcript: &mut Transcript,
    ) -> SumcheckOutput<F> {
        let num_rounds = Self::assert_shape(&polys, degree);

        let mut round_polys = Vec::with_capacity(num_rounds);
        let mut r = Vec::with_capacity(num_rounds);

        for _ in 0..num_rounds {
            let poly = Self::round_polynomial(&polys, claim, degree, &comb);
            let coeffs = poly.to_coeffs();

            transcript.append_serializable(b"round_poly", &coeffs);
            let r_i = transcript.challenge_field::<F>(b"challenge_r");
            claim = coeffs.eval(r_i);

            for p in polys.iter_mut() {
                p.fold(r_i);
            }

            round_polys.push(poly);
            r.push(r_i);
        }

        let final_evals = polys.iter().map(|p| p.final_constant()).collect();

        SumcheckOutput {
            proof: Self { polys: round_polys },
            point: r,
            final_evals,
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr as F;

    use crate::primitive::{DenseMultilinearPoly, Transcript};

    use super::{SumcheckOutput, SumcheckProof};

    #[test]
    fn arbitrary_quadratic_sumcheck_with_nonzero_claim() {
        // Spartan step 2 の形: sum_y (r_a * A(y) + r_b * B(y)) * z(y)、次数 2。
        let (r_a, r_b) = (F::from(31), F::from(37));
        let a_evals = [2u64, 3, 5, 7, 11, 13, 17, 19].map(F::from);
        let b_evals = [23u64, 29, 31, 37, 41, 43, 47, 53].map(F::from);
        let z_evals = [1u64, 4, 9, 16, 25, 36, 49, 64].map(F::from);

        let comb = |v: &[F]| (r_a * v[0] + r_b * v[1]) * v[2];
        let expected_sum: F = (0..8)
            .map(|i| comb(&[a_evals[i], b_evals[i], z_evals[i]]))
            .sum();
        assert_ne!(expected_sum, F::from(0));

        let polys = vec![
            DenseMultilinearPoly::new(a_evals.to_vec(), 3),
            DenseMultilinearPoly::new(b_evals.to_vec(), 3),
            DenseMultilinearPoly::new(z_evals.to_vec(), 3),
        ];

        let mut prover_transcript = Transcript::new(b"arbitrary-sumcheck-test");
        let SumcheckOutput {
            proof,
            point: prover_r,
            final_evals,
        } = SumcheckProof::prove(polys, expected_sum, 2, comb, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"arbitrary-sumcheck-test");
        let (final_claim, verifier_r) = proof
            .verify(expected_sum, 2, &mut verifier_transcript)
            .unwrap();

        assert_eq!(prover_r, verifier_r);
        assert_eq!(verifier_r.len(), 3);
        assert_eq!(final_evals.len(), 3);
        assert_eq!(final_claim, comb(&final_evals));
    }

    #[test]
    fn sparse_and_dense_polys_give_the_same_proof() {
        // comb(0, 0) != 0 になる形を混ぜて、全 zero ペアの補正項が効くことも確かめる。
        // ペア (2,3) と (4,5) は両因子とも zero。
        let a = DenseMultilinearPoly::new([0u64, 3, 0, 0, 0, 0, 0, 7].map(F::from).to_vec(), 3);
        let b = DenseMultilinearPoly::new([2u64, 0, 0, 0, 0, 0, 5, 0].map(F::from).to_vec(), 3);
        let comb = |v: &[F]| (v[0] + F::from(1)) * v[1] + F::from(11);
        let claim: F = (0..8).map(|i| comb(&[a.evals()[i], b.evals()[i]])).sum();

        let mut dense_transcript = Transcript::new(b"sparse-vs-dense");
        let dense = SumcheckProof::prove(
            vec![a.clone(), b.clone()],
            claim,
            2,
            comb,
            &mut dense_transcript,
        );

        let mut nnz_transcript = Transcript::new(b"sparse-vs-dense");
        let sparse = SumcheckProof::prove(
            vec![a.to_sparse(), b.to_sparse()],
            claim,
            2,
            comb,
            &mut nnz_transcript,
        );

        assert_eq!(dense.point, sparse.point);
        assert_eq!(dense.final_evals, sparse.final_evals);

        let mut verifier_transcript = Transcript::new(b"sparse-vs-dense");
        let (final_claim, _) = sparse
            .proof
            .verify(claim, 2, &mut verifier_transcript)
            .unwrap();
        assert_eq!(final_claim, comb(&sparse.final_evals));
    }

    #[test]
    fn wrong_claim_is_rejected() {
        let polys = vec![
            DenseMultilinearPoly::new([1u64, 2, 3, 4].map(F::from).to_vec(), 2),
            DenseMultilinearPoly::new([5u64, 6, 7, 8].map(F::from).to_vec(), 2),
        ];
        let comb = |v: &[F]| v[0] * v[1];
        let true_sum = F::from(5 + 12 + 21 + 32);

        let mut prover_transcript = Transcript::new(b"wrong-claim-test");
        let output = SumcheckProof::prove(polys, true_sum, 2, comb, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"wrong-claim-test");
        assert!(
            output
                .proof
                .verify(true_sum + F::from(1), 2, &mut verifier_transcript)
                .is_none()
        );
    }

    #[test]
    fn verifier_enforces_the_protocol_degree_bound() {
        let polys = vec![DenseMultilinearPoly::new(
            [1u64, 2, 3, 4].map(F::from).to_vec(),
            2,
        )];
        let claim = F::from(1 + 2 + 3 + 4);
        let mut prover_transcript = Transcript::new(b"degree-bound-test");
        // prover が次数 3 のラウンド多項式を送っても、プロトコルの上限が 2 なら拒否される。
        let output = SumcheckProof::prove(polys, claim, 3, |v| v[0], &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"degree-bound-test");
        assert!(
            output
                .proof
                .verify(claim, 2, &mut verifier_transcript)
                .is_none()
        );
    }
}
