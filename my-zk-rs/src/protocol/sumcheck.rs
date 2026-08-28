use ark_ff::PrimeField;

use crate::primitive::{DenseMultilinearPoly, EvalsUniPoly, MultilinearPoly, Transcript, UniPoly};

pub struct SumcheckProof<F: PrimeField> {
    polys: Vec<EvalsUniPoly<F>>,
    max_degree: usize,
}

impl<F: PrimeField> SumcheckProof<F> {
    pub fn verify(proof: &SumcheckProof<F>, sum: F, transcript: &mut Transcript) -> (F, Vec<F>) {
        let mut h = sum;
        let mut r: Vec<F> = Vec::with_capacity(proof.polys.len());

        for p in &proof.polys {
            assert!(p.degree() <= proof.max_degree);

            let p = &p.to_coeffs();
            let eval_zero = p.eval(F::zero());
            let eval_one = p.eval(F::one());
            assert_eq!(eval_one + eval_zero, h);

            transcript.append_serializable(b"round_poly", p);
            let r_i = transcript.challenge_field::<F>(b"challenge_r");
            r.push(r_i);

            h = p.eval(r_i);
        }

        (h, r)
    }

    /// `eq * (Az * Bz - Cz)` の、初期 claim を zero とする cubic sum-check。
    ///
    /// 戻り値の最後の要素は `[eq(r), Az(r), Bz(r), Cz(r)]`。
    pub fn prove_step1_sumcheck(
        mut eq: DenseMultilinearPoly<F>,
        mut az: DenseMultilinearPoly<F>,
        mut bz: DenseMultilinearPoly<F>,
        mut cz: DenseMultilinearPoly<F>,
        transcript: &mut Transcript,
    ) -> (Self, Vec<F>, Vec<F>) {
        let num_rounds = eq.vars();
        assert_eq!(az.vars(), num_rounds, "eq and Az dimensions differ");
        assert_eq!(bz.vars(), num_rounds, "eq and Bz dimensions differ");
        assert_eq!(cz.vars(), num_rounds, "eq and Cz dimensions differ");

        let mut claim = F::zero();
        let mut polys = Vec::with_capacity(num_rounds);
        let mut r = Vec::with_capacity(num_rounds);

        for _ in 0..num_rounds {
            let mut eval_zero = F::zero();
            let mut eval_two = F::zero();
            let mut eval_three = F::zero();

            for i in 0..eq.evals().len() / 2 {
                let lo = 2 * i;
                let hi = lo + 1;

                let eq_zero = eq.evals()[lo];
                let az_zero = az.evals()[lo];
                let bz_zero = bz.evals()[lo];
                let cz_zero = cz.evals()[lo];
                eval_zero += eq_zero * (az_zero * bz_zero - cz_zero);

                let eq_two = eq.evals()[hi] + eq.evals()[hi] - eq_zero;
                let az_two = az.evals()[hi] + az.evals()[hi] - az_zero;
                let bz_two = bz.evals()[hi] + bz.evals()[hi] - bz_zero;
                let cz_two = cz.evals()[hi] + cz.evals()[hi] - cz_zero;
                eval_two += eq_two * (az_two * bz_two - cz_two);

                let eq_three = eq_two + eq.evals()[hi] - eq_zero;
                let az_three = az_two + az.evals()[hi] - az_zero;
                let bz_three = bz_two + bz.evals()[hi] - bz_zero;
                let cz_three = cz_two + cz.evals()[hi] - cz_zero;
                eval_three += eq_three * (az_three * bz_three - cz_three);
            }

            let poly = EvalsUniPoly::new(vec![eval_zero, claim - eval_zero, eval_two, eval_three]);
            let coeffs = poly.to_coeffs();

            transcript.append_serializable(b"round_poly", &coeffs);
            let r_i = transcript.challenge_field::<F>(b"challenge_r");
            claim = coeffs.eval(r_i);

            eq.fold(r_i);
            az.fold(r_i);
            bz.fold(r_i);
            cz.fold(r_i);

            polys.push(poly);
            r.push(r_i);
        }

        let claims = vec![
            eq.final_constant(),
            az.final_constant(),
            bz.final_constant(),
            cz.final_constant(),
        ];

        (
            Self {
                polys,
                max_degree: 3,
            },
            r,
            claims,
        )
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr as F;

    use crate::primitive::{DenseMultilinearPoly, EqPoly, Transcript};

    use super::SumcheckProof;

    #[test]
    fn step1_prover_and_verifier_agree() {
        let eq = EqPoly::new(vec![F::from(23), F::from(29)]).to_dense_mlp();
        let az_evals = [2u64, 3, 5, 7].map(F::from);
        let bz_evals = [11u64, 13, 17, 19].map(F::from);
        let cz_evals = az_evals
            .iter()
            .zip(&bz_evals)
            .map(|(&az, &bz)| az * bz)
            .collect();

        let az = DenseMultilinearPoly::new(az_evals.to_vec(), 2);
        let bz = DenseMultilinearPoly::new(bz_evals.to_vec(), 2);
        let cz = DenseMultilinearPoly::new(cz_evals, 2);

        let mut prover_transcript = Transcript::new(b"step1-sumcheck-test");
        let (proof, prover_r, claims) =
            SumcheckProof::prove_step1_sumcheck(eq, az, bz, cz, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"step1-sumcheck-test");
        let (final_claim, verifier_r) =
            SumcheckProof::verify(&proof, F::from(0), &mut verifier_transcript);

        assert_eq!(prover_r, verifier_r);
        assert_eq!(claims.len(), 4);
        assert_eq!(final_claim, claims[0] * (claims[1] * claims[2] - claims[3]));
    }
}
