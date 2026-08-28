use std::marker::PhantomData;

use ark_ff::Field;

use crate::primitive::Matrix;

/// R1CS の変数配置。assignment は `z = [io..., 1, witness..., 0...]`。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct R1CSStructure<F: Field> {
    pub num_constraints: usize,
    pub num_io: usize,
    pub num_witness: usize,
    _f: PhantomData<F>,
}

impl<F: Field> R1CSStructure<F> {
    pub fn new(num_constraints: usize, num_io: usize, num_witness: usize) -> Self {
        Self {
            num_constraints,
            num_io,
            num_witness,
            _f: PhantomData,
        }
    }
}

/// 制約系 $A z \circ B z = C z$。
pub struct R1CS<F: Field, M: Matrix<F>> {
    pub a: M,
    pub b: M,
    pub c: M,
    pub structure: R1CSStructure<F>,
}

impl<F: Field, M: Matrix<F>> R1CS<F, M> {
    pub fn new(a: M, b: M, c: M, structure: R1CSStructure<F>) -> Self {
        assert!(
            a.shape() == b.shape() && a.shape() == c.shape(),
            "R1CS matrices must share the same shape"
        );
        assert!(
            a.rows().is_power_of_two() && a.cols().is_power_of_two(),
            "R1CS matrices must have power-of-two dimensions"
        );
        assert!(
            structure.num_constraints <= a.rows(),
            "too many constraints for the matrices"
        );
        assert!(
            structure.num_vars <= a.cols(),
            "too many variables for the matrices"
        );

        Self { a, b, c, structure }
    }

    /// `z = [io..., 1, witness..., 0...]` を行列の列数まで zero-padding して作る。
    pub fn assignment(&self, io: &[F], witness: &[F]) -> Vec<F> {
        assert_eq!(io.len(), self.structure.num_io, "invalid io length");
        assert_eq!(
            witness.len(),
            self.structure.num_witness(),
            "invalid witness length"
        );

        let one_index = self.structure.num_io;
        let witness_start = one_index + 1;

        let mut z = vec![F::zero(); self.a.cols()];
        z[..io.len()].copy_from_slice(io);
        z[one_index] = F::one();
        z[witness_start..self.structure.num_vars].copy_from_slice(witness);
        z
    }

    /// $A z \circ B z = C z$ が最初の `num_constraints` 行で成り立つか。
    pub fn is_sat(&self, z: &[F]) -> bool {
        assert_eq!(z.len(), self.a.cols(), "invalid assignment length");

        let az = self.a.mul_vec(z);
        let bz = self.b.mul_vec(z);
        let cz = self.c.mul_vec(z);

        az.iter()
            .zip(&bz)
            .zip(&cz)
            .take(self.structure.num_constraints)
            .all(|((&a, &b), &c)| a * b == c)
    }
}

#[cfg(test)]
mod tests {
    use super::{R1CS, R1CSStructure};
    use crate::primitive::DenseMatrix;
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};

    #[test]
    fn assignment_and_satisfaction_use_the_r1cs_layout() {
        type F = Fr;

        // z = [x, 1, y, out], with x * y = out. The matrices have
        // four extra zero columns to model the MLE-friendly padding.
        let r1cs = R1CS::new(
            DenseMatrix::from_usize([[1, 0, 0, 0, 0, 0, 0, 0]]),
            DenseMatrix::from_usize([[0, 0, 1, 0, 0, 0, 0, 0]]),
            DenseMatrix::from_usize([[0, 0, 0, 1, 0, 0, 0, 0]]),
            R1CSStructure::new(1, 4, 1),
        );

        let assignment = r1cs.assignment(&[F::from(3u64)], &[F::from(5u64), F::from(15u64)]);
        assert_eq!(
            assignment,
            [
                F::from(3u64),
                F::one(),
                F::from(5u64),
                F::from(15u64),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
            ]
        );
        assert!(r1cs.is_sat(&assignment));

        let unsatisfied = r1cs.assignment(&[F::from(3u64)], &[F::from(5u64), F::from(14u64)]);
        assert!(!r1cs.is_sat(&unsatisfied));
    }
}
