#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

use std::{hint::black_box, time::Duration};

use ark_ff::Zero;
use ark_poly::{
    DenseMVPolynomial, Polynomial,
    polynomial::multivariate::{SparsePolynomial, SparseTerm, Term},
};
use criterion::{Criterion, criterion_group, criterion_main};
use my_zk_rs::{
    Fp25519,
    primitive::{Monomial, MvPolynomial},
};

const VARS: usize = 4;
const TERMS: usize = 128;
const MAX_EXP: usize = 5;

type MyPoly = MvPolynomial<Fp25519, VARS>;
type ArkPoly = SparsePolynomial<Fp25519, SparseTerm>;

fn f(x: u64) -> Fp25519 {
    Fp25519::from_u64(x)
}

fn next_u64(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(0xda94_2042_e4dd_58b5)
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    *state
}

fn sample_exponents(seed: u64, term_index: usize) -> [usize; VARS] {
    let mut state = seed ^ ((term_index as u64) << 32);
    let mut exps = [0usize; VARS];

    for exp in &mut exps {
        *exp = (next_u64(&mut state) as usize) % (MAX_EXP + 1);
    }

    if exps.iter().all(|&exp| exp == 0) {
        exps[term_index % VARS] = 1;
    }

    exps
}

fn sample_coeff(seed: u64, term_index: usize) -> Fp25519 {
    f(seed.wrapping_add((term_index as u64 + 1) * 17))
}

fn ark_term(exps: [usize; VARS]) -> SparseTerm {
    SparseTerm::new(
        exps.into_iter()
            .enumerate()
            .filter_map(|(var, exp)| (exp != 0).then_some((var, exp)))
            .collect(),
    )
}

fn ark_term_mul(lhs: &SparseTerm, rhs: &SparseTerm) -> SparseTerm {
    SparseTerm::new(lhs.iter().chain(rhs.iter()).copied().collect())
}

fn ark_mul(lhs: &ArkPoly, rhs: &ArkPoly) -> ArkPoly {
    if lhs.is_zero() || rhs.is_zero() {
        return ArkPoly::zero();
    }

    let mut terms = Vec::with_capacity(lhs.terms.len().saturating_mul(rhs.terms.len()));
    for (lhs_coeff, lhs_term) in &lhs.terms {
        for (rhs_coeff, rhs_term) in &rhs.terms {
            terms.push((*lhs_coeff * *rhs_coeff, ark_term_mul(lhs_term, rhs_term)));
        }
    }

    ArkPoly::from_coefficients_vec(lhs.num_vars.max(rhs.num_vars), terms)
}

fn make_polys(seed: u64) -> (MyPoly, ArkPoly) {
    let my_terms = make_my_terms(seed);
    let ark_terms = make_ark_terms(&my_terms);

    (
        MyPoly::from_terms(my_terms),
        ArkPoly::from_coefficients_vec(VARS, ark_terms),
    )
}

fn make_my_terms(seed: u64) -> Vec<(Monomial<VARS>, Fp25519)> {
    (0..TERMS)
        .map(|i| {
            let exps = sample_exponents(seed, i);
            (Monomial::new(exps), sample_coeff(seed, i))
        })
        .collect()
}

fn make_ark_terms(my_terms: &[(Monomial<VARS>, Fp25519)]) -> Vec<(Fp25519, SparseTerm)> {
    my_terms
        .iter()
        .map(|(monomial, coeff)| (*coeff, ark_term(*monomial.exponents())))
        .collect()
}

fn sample_point() -> [Fp25519; VARS] {
    [f(3), f(5), f(7), f(11)]
}

fn assert_same_behavior(my_lhs: &MyPoly, my_rhs: &MyPoly, ark_lhs: &ArkPoly, ark_rhs: &ArkPoly) {
    let point = sample_point();
    let ark_point = point.to_vec();

    assert_eq!(my_lhs.degree().unwrap_or_default(), ark_lhs.degree());
    assert_eq!(my_rhs.degree().unwrap_or_default(), ark_rhs.degree());
    assert_eq!(my_lhs.eval(&point), ark_lhs.evaluate(&ark_point));
    assert_eq!(my_rhs.eval(&point), ark_rhs.evaluate(&ark_point));
    assert_eq!(
        (my_lhs + my_rhs).eval(&point),
        (ark_lhs + ark_rhs).evaluate(&ark_point)
    );
    assert_eq!(
        (my_lhs - my_rhs).eval(&point),
        (ark_lhs - ark_rhs).evaluate(&ark_point)
    );
    assert_eq!(
        (my_lhs * my_rhs).eval(&point),
        ark_mul(ark_lhs, ark_rhs).evaluate(&ark_point)
    );
}

fn bench_poly_vs_ark(c: &mut Criterion) {
    let (my_lhs, ark_lhs) = make_polys(0x243f_6a88_85a3_08d3);
    let (my_rhs, ark_rhs) = make_polys(0x1319_8a2e_0370_7344);
    assert_same_behavior(&my_lhs, &my_rhs, &ark_lhs, &ark_rhs);

    let init_my_terms = make_my_terms(0xa409_3822_299f_31d0);
    let init_ark_terms = make_ark_terms(&init_my_terms);
    let point = sample_point();
    let ark_point = point.to_vec();
    let label = format!("vars{VARS}_terms{TERMS}_max_exp{MAX_EXP}");

    let mut group = c.benchmark_group(format!("poly_init/{label}"));
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_zk_rs", |b| {
        b.iter(|| MyPoly::from_terms(black_box(init_my_terms.clone())))
    });
    group.bench_function("ark_poly", |b| {
        b.iter(|| ArkPoly::from_coefficients_vec(VARS, black_box(init_ark_terms.clone())))
    });
    group.finish();

    let mut group = c.benchmark_group(format!("poly_eval/{label}"));
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_zk_rs", |b| {
        b.iter(|| black_box(&my_lhs).eval(black_box(&point)))
    });
    group.bench_function("ark_poly", |b| {
        b.iter(|| black_box(&ark_lhs).evaluate(black_box(&ark_point)))
    });
    group.finish();

    let mut group = c.benchmark_group(format!("poly_add/{label}"));
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_zk_rs", |b| {
        b.iter(|| black_box(black_box(&my_lhs) + black_box(&my_rhs)))
    });
    group.bench_function("ark_poly", |b| {
        b.iter(|| black_box(black_box(&ark_lhs) + black_box(&ark_rhs)))
    });
    group.finish();

    let mut group = c.benchmark_group(format!("poly_sub/{label}"));
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_zk_rs", |b| {
        b.iter(|| black_box(black_box(&my_lhs) - black_box(&my_rhs)))
    });
    group.bench_function("ark_poly", |b| {
        b.iter(|| black_box(black_box(&ark_lhs) - black_box(&ark_rhs)))
    });
    group.finish();

    let mut group = c.benchmark_group(format!("poly_mul/{label}"));
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_zk_rs", |b| {
        b.iter(|| black_box(black_box(&my_lhs) * black_box(&my_rhs)))
    });
    group.bench_function("ark_poly_naive", |b| {
        b.iter(|| black_box(ark_mul(black_box(&ark_lhs), black_box(&ark_rhs))))
    });
    group.finish();

    let mut group = c.benchmark_group(format!("poly_degree/{label}"));
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_zk_rs", |b| b.iter(|| black_box(&my_lhs).degree()));
    group.bench_function("ark_poly", |b| b.iter(|| black_box(&ark_lhs).degree()));
    group.finish();
}

criterion_group!(benches, bench_poly_vs_ark);
criterion_main!(benches);
