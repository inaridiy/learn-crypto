#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use my_zk_rs::{BoolHyperCube, Fp25519, MvPolynomial, multilinear_extension};
use sha2::{Digest, Sha256};
use std::fmt;

struct FiatShamirTranscript {
    state: Sha256,
}

impl FiatShamirTranscript {
    fn new(domain_separator: &'static [u8]) -> Self {
        let mut transcript = Self {
            state: Sha256::new(),
        };
        transcript.append_bytes(b"domain", domain_separator);
        transcript
    }

    fn append_usize(&mut self, label: &'static [u8], value: usize) {
        self.append_bytes(label, &value.to_le_bytes());
    }

    fn append_bytes(&mut self, label: &'static [u8], bytes: &[u8]) {
        self.state.update((label.len() as u64).to_le_bytes());
        self.state.update(label);
        self.state.update((bytes.len() as u64).to_le_bytes());
        self.state.update(bytes);
    }

    fn append_field(&mut self, label: &'static [u8], value: &Fp25519) {
        let mut bytes = Vec::with_capacity(value.uncompressed_size());
        value
            .serialize_uncompressed(&mut bytes)
            .expect("serializing a field element into Vec cannot fail");
        self.append_bytes(label, &bytes);
    }

    fn append_poly<const N: usize>(
        &mut self,
        label: &'static [u8],
        poly: &MvPolynomial<Fp25519, N>,
    ) {
        self.append_bytes(b"poly-label", label);
        self.append_usize(b"poly-num-vars", N);
        self.append_usize(b"poly-num-terms", poly.num_terms());

        for (monomial, coeff) in poly.terms() {
            for exp in monomial.exponents() {
                self.append_usize(b"monomial-exp", *exp);
            }
            self.append_field(b"coeff", coeff);
        }
    }

    fn append_univariate_poly(&mut self, label: &'static [u8], poly: &UnivariatePolynomial) {
        self.append_bytes(b"univariate-poly-label", label);
        self.append_usize(b"univariate-poly-degree-bound", poly.coeffs.len() - 1);
        for coeff in &poly.coeffs {
            self.append_field(b"univariate-poly-coeff", coeff);
        }
    }

    fn challenge_field(&mut self, label: &'static [u8]) -> Fp25519 {
        for counter in 0u64.. {
            let mut hasher = self.state.clone();
            hasher.update((label.len() as u64).to_le_bytes());
            hasher.update(label);
            hasher.update(counter.to_le_bytes());
            let digest = hasher.finalize();

            if let Some(challenge) = Fp25519::from_random_bytes(&digest) {
                self.append_field(label, &challenge);
                return challenge;
            }
        }

        unreachable!("unbounded counter eventually samples a field element")
    }
}

#[derive(Clone, Debug)]
struct UnivariatePolynomial {
    coeffs: Vec<Fp25519>,
}

impl UnivariatePolynomial {
    fn quadratic_from_evaluations_at_0_1_2(y0: Fp25519, y1: Fp25519, y2: Fp25519) -> Self {
        let two = Fp25519::from(2);
        let c2 = (y2 - two * y1 + y0) / two;
        let c1 = y1 - y0 - c2;

        Self {
            coeffs: vec![y0, c1, c2],
        }
    }

    fn eval(&self, x: Fp25519) -> Fp25519 {
        self.coeffs
            .iter()
            .rev()
            .fold(Fp25519::zero(), |value, coeff| value * x + coeff)
    }

    fn sum_over_boolean_hypercube(&self) -> Fp25519 {
        self.eval(Fp25519::zero()) + self.eval(Fp25519::one())
    }
}

impl fmt::Display for UnivariatePolynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut wrote_term = false;

        for (degree, coeff) in self.coeffs.iter().enumerate() {
            if *coeff == Fp25519::zero() {
                continue;
            }

            if wrote_term {
                write!(f, " + ")?;
            }

            match degree {
                0 => write!(f, "{coeff}")?,
                1 if *coeff == Fp25519::one() => write!(f, "t")?,
                1 => write!(f, "{coeff}*t")?,
                _ if *coeff == Fp25519::one() => write!(f, "t^{degree}")?,
                _ => write!(f, "{coeff}*t^{degree}")?,
            }

            wrote_term = true;
        }

        if wrote_term { Ok(()) } else { write!(f, "0") }
    }
}

#[derive(Debug)]
struct NonInteractiveSumCheckProof {
    claimed_sum: Fp25519,
    round_polys: Vec<UnivariatePolynomial>,
    challenges: Vec<Fp25519>,
    final_evaluation: Fp25519,
}

fn append_sumcheck_statement(
    transcript: &mut FiatShamirTranscript,
    a: &[Fp25519; 5],
    b: &[Fp25519; 5],
    c: &[Fp25519; 5],
    g: &MvPolynomial<Fp25519, 3>,
    claimed_sum: Fp25519,
) {
    for value in a {
        transcript.append_field(b"a", value);
    }
    for value in b {
        transcript.append_field(b"b", value);
    }
    for value in c {
        transcript.append_field(b"c", value);
    }
    transcript.append_poly(b"g", g);
    transcript.append_field(b"claimed-sum", &claimed_sum);
}

fn prove_sumcheck(
    g: &MvPolynomial<Fp25519, 3>,
    claimed_sum: Fp25519,
    transcript: &mut FiatShamirTranscript,
) -> NonInteractiveSumCheckProof {
    let mut round_polys = Vec::with_capacity(3);
    let mut challenges = Vec::with_capacity(3);

    for round in 0..3 {
        let round_poly = sumcheck_round_poly(g, &challenges, round);
        transcript.append_usize(b"sumcheck-round", round);
        transcript.append_univariate_poly(b"sumcheck-round-poly", &round_poly);
        let challenge = transcript.challenge_field(b"sumcheck-challenge");

        round_polys.push(round_poly);
        challenges.push(challenge);
    }

    let final_point = [challenges[0], challenges[1], challenges[2]];
    let final_evaluation = g.eval(&final_point);

    NonInteractiveSumCheckProof {
        claimed_sum,
        round_polys,
        challenges,
        final_evaluation,
    }
}

fn verify_sumcheck(
    g: &MvPolynomial<Fp25519, 3>,
    expected_claim: Fp25519,
    proof: &NonInteractiveSumCheckProof,
    transcript: &mut FiatShamirTranscript,
) -> bool {
    if proof.claimed_sum != expected_claim
        || proof.round_polys.len() != 3
        || proof.challenges.len() != 3
    {
        return false;
    }

    let mut expected_sum = proof.claimed_sum;
    let mut challenges = Vec::with_capacity(3);

    for (round, round_poly) in proof.round_polys.iter().enumerate() {
        if round_poly.sum_over_boolean_hypercube() != expected_sum {
            return false;
        }

        transcript.append_usize(b"sumcheck-round", round);
        transcript.append_univariate_poly(b"sumcheck-round-poly", round_poly);
        let challenge = transcript.challenge_field(b"sumcheck-challenge");

        if proof.challenges[round] != challenge {
            return false;
        }

        expected_sum = round_poly.eval(challenge);
        challenges.push(challenge);
    }

    let final_point = [challenges[0], challenges[1], challenges[2]];
    proof.final_evaluation == expected_sum && g.eval(&final_point) == proof.final_evaluation
}

fn sumcheck_round_poly(
    g: &MvPolynomial<Fp25519, 3>,
    fixed_prefix: &[Fp25519],
    round: usize,
) -> UnivariatePolynomial {
    let y0 = sum_with_round_variable(g, fixed_prefix, round, Fp25519::zero());
    let y1 = sum_with_round_variable(g, fixed_prefix, round, Fp25519::one());
    let y2 = sum_with_round_variable(g, fixed_prefix, round, Fp25519::from(2));

    UnivariatePolynomial::quadratic_from_evaluations_at_0_1_2(y0, y1, y2)
}

fn sum_with_round_variable(
    g: &MvPolynomial<Fp25519, 3>,
    fixed_prefix: &[Fp25519],
    round: usize,
    round_value: Fp25519,
) -> Fp25519 {
    let remaining_vars = 3 - round - 1;
    let mut sum = Fp25519::zero();

    for mask in 0..(1 << remaining_vars) {
        let mut point = [Fp25519::zero(); 3];

        point[..round].copy_from_slice(&fixed_prefix[..round]);
        point[round] = round_value;

        for i in 0..remaining_vars {
            point[round + 1 + i] = if ((mask >> i) & 1) == 1 {
                Fp25519::one()
            } else {
                Fp25519::zero()
            };
        }

        sum += g.eval(&point);
    }

    sum
}

fn sum_over_boolean_hypercube(g: &MvPolynomial<Fp25519, 3>) -> Fp25519 {
    BoolHyperCube::<3>::iter()
        .map(|x| g.eval(&x.to_field_point()))
        .sum()
}

fn main() {
    type F = Fp25519;

    let a = [10, 20, 30, 40, 50].map(F::from);
    let b = [10, 20, 30, 40, 50].map(F::from);
    let c = [100, 400, 900, 1600, 2500].map(F::from);

    let ta = multilinear_extension(a);
    let tb = multilinear_extension(b);
    let tc = multilinear_extension(c);

    println!("a MLE: {ta}");
    println!("b MLE: {tb}");
    println!("c MLE: {tc}");

    // g(x) := \tilde{a}(x) \dot  \tilde{b}(x) -  \tilde{c}(x)
    let g = ta * tb - tc;

    let mut q = MvPolynomial::<F, 3>::zero();
    for x in BoolHyperCube::<3>::iter() {
        q += x.teq().scale(g.eval(&x.to_field_point()))
    }

    let mut transcript = FiatShamirTranscript::new(b"a_dot_b_eq_c/v1");
    for value in &a {
        transcript.append_field(b"a", value);
    }
    for value in &b {
        transcript.append_field(b"b", value);
    }
    for value in &c {
        transcript.append_field(b"c", value);
    }
    transcript.append_poly(b"q", &q);

    let r = [
        transcript.challenge_field(b"r[0]"),
        transcript.challenge_field(b"r[1]"),
        transcript.challenge_field(b"r[2]"),
    ];
    let h = q.eval(&r);

    println!("Fiat-Shamir challenge r: [{}, {}, {}]", r[0], r[1], r[2]);
    println!("h {h}");

    let claimed_sum = sum_over_boolean_hypercube(&g);
    let mut prover_transcript = FiatShamirTranscript::new(b"a_dot_b_eq_c/sumcheck/v1");
    append_sumcheck_statement(&mut prover_transcript, &a, &b, &c, &g, claimed_sum);
    let sumcheck_proof = prove_sumcheck(&g, claimed_sum, &mut prover_transcript);

    let mut verifier_transcript = FiatShamirTranscript::new(b"a_dot_b_eq_c/sumcheck/v1");
    append_sumcheck_statement(&mut verifier_transcript, &a, &b, &c, &g, claimed_sum);
    let sumcheck_verified =
        verify_sumcheck(&g, claimed_sum, &sumcheck_proof, &mut verifier_transcript);

    println!("SumCheck claimed sum: {}", sumcheck_proof.claimed_sum);
    for (round, round_poly) in sumcheck_proof.round_polys.iter().enumerate() {
        println!("SumCheck round {} polynomial: {}", round + 1, round_poly);
    }
    println!(
        "SumCheck challenges: [{}, {}, {}]",
        sumcheck_proof.challenges[0], sumcheck_proof.challenges[1], sumcheck_proof.challenges[2]
    );
    println!(
        "SumCheck final evaluation g(r): {}",
        sumcheck_proof.final_evaluation
    );
    println!("SumCheck verified: {sumcheck_verified}");
}
