use ark_ff::Field;

use super::helpers::log2_ceil;
use super::matrix::Matrix;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct R1CSStructure {
    pub num_constraints: usize,
    pub num_vars: usize,
    pub num_io: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct R1CS<F: Field, const CONSTRAINT_BITS: usize = 0, const VAR_BITS: usize = 0>
where
    [(); 1 << CONSTRAINT_BITS]:,
    [(); 1 << VAR_BITS]:,
{
    pub a: Matrix<F, CONSTRAINT_BITS, VAR_BITS>,
    pub b: Matrix<F, CONSTRAINT_BITS, VAR_BITS>,
    pub c: Matrix<F, CONSTRAINT_BITS, VAR_BITS>,
    pub structure: R1CSStructure,
}

impl<F: Field, const CONSTRAINT_BITS: usize, const VAR_BITS: usize>
    R1CS<F, CONSTRAINT_BITS, VAR_BITS>
where
    [(); 1 << CONSTRAINT_BITS]:,
    [(); 1 << VAR_BITS]:,
{
    pub fn new(
        a: Matrix<F, CONSTRAINT_BITS, VAR_BITS>,
        b: Matrix<F, CONSTRAINT_BITS, VAR_BITS>,
        c: Matrix<F, CONSTRAINT_BITS, VAR_BITS>,
        structure: R1CSStructure,
    ) -> Self {
        assert!(
            structure.num_constraints <= 1 << CONSTRAINT_BITS,
            "too many constraints"
        );
        assert!(structure.num_vars <= 1 << VAR_BITS, "too many variables");
        assert!(
            structure.num_io < structure.num_vars,
            "R1CS structure does not fit"
        );
        Self { a, b, c, structure }
    }

    pub fn from_unpadded<const NUM_CONSTRAINTS: usize, const NUM_VARS: usize>(
        a: [[F; NUM_VARS]; NUM_CONSTRAINTS],
        b: [[F; NUM_VARS]; NUM_CONSTRAINTS],
        c: [[F; NUM_VARS]; NUM_CONSTRAINTS],
        num_io: usize,
    ) -> R1CS<F, { log2_ceil(NUM_CONSTRAINTS) }, { log2_ceil(NUM_VARS) }>
    where
        [(); 1 << log2_ceil(NUM_CONSTRAINTS)]:,
        [(); 1 << log2_ceil(NUM_VARS)]:,
    {
        R1CS::<F, { log2_ceil(NUM_CONSTRAINTS) }, { log2_ceil(NUM_VARS) }>::new(
            Matrix::<F>::from_unpadded(a),
            Matrix::<F>::from_unpadded(b),
            Matrix::<F>::from_unpadded(c),
            R1CSStructure {
                num_constraints: NUM_CONSTRAINTS,
                num_vars: NUM_VARS,
                num_io,
            },
        )
    }

    pub fn from_unpadded_usize<const NUM_CONSTRAINTS: usize, const NUM_VARS: usize>(
        a: [[usize; NUM_VARS]; NUM_CONSTRAINTS],
        b: [[usize; NUM_VARS]; NUM_CONSTRAINTS],
        c: [[usize; NUM_VARS]; NUM_CONSTRAINTS],
        num_io: usize,
    ) -> R1CS<F, { log2_ceil(NUM_CONSTRAINTS) }, { log2_ceil(NUM_VARS) }>
    where
        [(); 1 << log2_ceil(NUM_CONSTRAINTS)]:,
        [(); 1 << log2_ceil(NUM_VARS)]:,
    {
        Self::from_unpadded(
            a.map(|row| row.map(|value| F::from(value as u128))),
            b.map(|row| row.map(|value| F::from(value as u128))),
            c.map(|row| row.map(|value| F::from(value as u128))),
            num_io,
        )
    }

    pub fn assignment(&self, io: &[F], witness: &[F]) -> [F; 1 << VAR_BITS] {
        let one_index = self.structure.num_io;
        let witness_start = one_index + 1;

        assert_eq!(io.len(), self.structure.num_io, "invalid io length");
        assert_eq!(
            witness.len(),
            self.structure.num_vars - witness_start,
            "invalid witness length"
        );

        let mut assignment = [F::zero(); 1 << VAR_BITS];
        assignment[..io.len()].copy_from_slice(io);
        assignment[one_index] = F::one();
        assignment[witness_start..self.structure.num_vars].copy_from_slice(witness);
        assignment
    }

    pub fn is_sat(&self, assignment: &[F; 1 << VAR_BITS]) -> bool {
        (0..self.structure.num_constraints).all(|constraint_index| {
            self.a.dot_row(constraint_index, assignment)
                * self.b.dot_row(constraint_index, assignment)
                == self.c.dot_row(constraint_index, assignment)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{Matrix, R1CS};
    use ark_bls12_381::Fr as F;

    fn f(value: u64) -> F {
        F::from(value)
    }

    #[test]
    fn matrix_multiplies_by_vector() {
        let matrix = Matrix::<F>::from_unpadded([[f(1), f(2), f(3)], [f(4), f(5), f(6)]]);

        assert_eq!(matrix.rows().len(), 2);
        assert_eq!(matrix.row(0), &[f(1), f(2), f(3), f(0)]);
        assert_eq!(
            matrix.mul_vector(&[f(7), f(11), f(13), f(0)]),
            [f(68), f(161)]
        );
    }

    #[test]
    fn r1cs_accepts_satisfying_witness() {
        // z = [1, x, y, output], constraint: x * y = output.
        let r1cs = R1CS::<F>::from_unpadded(
            [[f(0), f(1), f(0), f(0)]],
            [[f(0), f(0), f(1), f(0)]],
            [[f(0), f(0), f(0), f(1)]],
            0,
        );

        assert!(r1cs.is_sat(&[f(1), f(3), f(5), f(15)]));
    }

    #[test]
    fn r1cs_rejects_unsatisfying_witness() {
        let r1cs = R1CS::<F>::from_unpadded(
            [[f(0), f(1), f(0), f(0)]],
            [[f(0), f(0), f(1), f(0)]],
            [[f(0), f(0), f(0), f(1)]],
            0,
        );

        assert!(!r1cs.is_sat(&[f(1), f(3), f(5), f(14)]));
    }
}
