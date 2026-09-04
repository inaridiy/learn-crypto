use std::marker::PhantomData;

use ark_ff::Field;

use crate::primitive::Matrix;

/// R1CS の変数個数。`z` の並びは制約系の型([`R1CS`] など)が定める。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct R1CSStructure {
    pub num_constraints: usize,
    pub num_io: usize,
    pub num_witness: usize,
}

impl R1CSStructure {
    pub fn new(num_constraints: usize, num_io: usize, num_witness: usize) -> Self {
        Self {
            num_constraints,
            num_io,
            num_witness,
        }
    }

    /// Pinocchio 流の `z = [1, io..., witness...]` の全長。
    pub fn dense_len(&self) -> usize {
        1 + self.num_io + self.num_witness
    }
}

/// 制約系 $A z \circ B z = C z$。Pinocchio 流の dense な配置
/// `z = [1, io..., witness...]` で、padding は行わない。
pub struct R1CS<F: Field, M: Matrix<F>> {
    pub a: M,
    pub b: M,
    pub c: M,
    pub structure: R1CSStructure,
    _f: PhantomData<F>,
}

impl<F: Field, M: Matrix<F>> R1CS<F, M> {
    pub fn new(a: M, b: M, c: M, structure: R1CSStructure) -> Self {
        assert!(
            a.shape() == b.shape() && a.shape() == c.shape(),
            "R1CS matrices must share the same shape"
        );
        assert_eq!(
            a.rows(),
            structure.num_constraints,
            "matrix rows must match the number of constraints"
        );
        assert_eq!(
            a.cols(),
            structure.dense_len(),
            "matrix columns must match the assignment length"
        );

        Self {
            a,
            b,
            c,
            structure,
            _f: PhantomData,
        }
    }

    /// `z = [1, io..., witness...]` を作る。
    pub fn assignment(&self, io: &[F], witness: &[F]) -> Vec<F> {
        assert_eq!(io.len(), self.structure.num_io, "invalid io length");
        assert_eq!(
            witness.len(),
            self.structure.num_witness,
            "invalid witness length"
        );

        let mut z = Vec::with_capacity(self.structure.dense_len());
        z.push(F::one());
        z.extend_from_slice(io);
        z.extend_from_slice(witness);
        z
    }

    /// $A z \circ B z = C z$ が全行で成り立つか。
    pub fn is_sat(&self, z: &[F]) -> bool {
        assert_eq!(z.len(), self.a.cols(), "invalid assignment length");

        let az = self.a.mul_vec(z);
        let bz = self.b.mul_vec(z);
        let cz = self.c.mul_vec(z);

        az.iter().zip(&bz).zip(&cz).all(|((&a, &b), &c)| a * b == c)
    }
}

#[cfg(test)]
mod tests {
    use super::{R1CSStructure, R1CS};
    use crate::primitive::DenseMatrix;
    use ark_bls12_381::Fr as F;
    use ark_ff::One;

    /// x * y = out。io = [x, out], witness = [y] で z = [1, x, out, y]。
    fn mul_r1cs() -> R1CS<F, DenseMatrix<F, 1, 4>> {
        R1CS::new(
            DenseMatrix::from_usize([[0, 1, 0, 0]]),
            DenseMatrix::from_usize([[0, 0, 0, 1]]),
            DenseMatrix::from_usize([[0, 0, 1, 0]]),
            R1CSStructure::new(1, 2, 1),
        )
    }

    #[test]
    fn dense_assignment_and_satisfaction() {
        let r1cs = mul_r1cs();

        let z = r1cs.assignment(&[F::from(3), F::from(15)], &[F::from(5)]);
        assert_eq!(z, [F::one(), F::from(3), F::from(15), F::from(5)]);
        assert!(r1cs.is_sat(&z));

        let bad = r1cs.assignment(&[F::from(3), F::from(14)], &[F::from(5)]);
        assert!(!r1cs.is_sat(&bad));
    }
}
