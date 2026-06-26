#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

use ark_ff::{Field, PrimeField};
use ark_serialize::CanonicalSerialize;
use my_zk_rs::primitive::{
    BoolHyperCube, Matrix, MvPolynomial, R1CS, Transcript, mle_from_hypercube_evaluations,
    mle_from_matrix,
};

fn init_transcript<F, const N: usize, const M: usize>(r1cs: &R1CS<F, N, M>) -> Transcript
where
    F: Field + CanonicalSerialize,
    [(); 1 << N]:,
    [(); 1 << M]:,
{
    let mut transcript = Transcript::new(b"spartan/nizk");
    transcript.append_bytes(b"field", std::any::type_name::<F>().as_bytes());
    transcript.append_usize(b"constraint-bits", N);
    transcript.append_usize(b"var-bits", M);
    transcript.append_usize(b"num-constraints", r1cs.structure.num_constraints);
    transcript.append_usize(b"num-vars", r1cs.structure.num_vars);
    transcript.append_usize(b"num-io", r1cs.structure.num_io);
    transcript.append_matrix(b"r1cs-a", &r1cs.a);
    transcript.append_matrix(b"r1cs-b", &r1cs.b);
    transcript.append_matrix(b"r1cs-c", &r1cs.c);
    transcript
}

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

pub fn prove<F, const N: usize, const M: usize>(r1cs: R1CS<F, N, M>, io: &[F], witness: &[F])
where
    F: CanonicalSerialize + PrimeField,
    [(); 1 << N]:,
    [(); 1 << M]:,
    [(); 1 << (N + M)]:,
{
    let assignment = r1cs.assignment(io, witness);

    assert!(r1cs.is_sat(&assignment), "R1CSがsatされません。");

    let mut transcript = init_transcript(&r1cs);

    let _rs = transcript.challenge_field::<F>(b"first-sumcheck");
    let tassignment = mle_from_hypercube_evaluations(assignment);
    let _bar_a = calc_bar_matrix(&r1cs.a, &tassignment);
}

pub fn main() {
    use ark_bls12_381::Fr as F;

    // out = x^3 + x + 5
    // w = [x,out,1,i1,i2,i3]

    let r1cs = R1CS::<F>::from_unpadded_usize(
        [
            [1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [1, 0, 0, 0, 1, 0],
            [0, 0, 5, 0, 0, 1],
        ],
        [
            [1, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
        ],
        [
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 1, 0, 0, 0, 0],
        ],
        2,
    );
    let assignment = r1cs.assignment(
        &[3, 35].map(F::from), //
        &[9, 27, 30].map(F::from),
    );

    assert!(r1cs.is_sat(&assignment));
}
