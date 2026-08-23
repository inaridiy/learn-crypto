use ark_ff::Field;

use super::matrix::Matrix;

/// assignment ベクトル $z$ (長さ $2^{\mathrm{VAR\_BITS}}$) は Spartan 論文に倣い、
/// 前半に public part、後半に witness を配置する:
///
/// ```text
/// z = ( io, 1, 0, ...  |  w, 0, ... )
///       前半 2^{VAR_BITS-1}   後半 2^{VAR_BITS-1}
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct R1CSStructure {
    pub num_constraints: usize,
    pub num_io: usize,
    pub num_witness: usize,
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
        assert!(VAR_BITS >= 1, "the assignment layout needs both halves");

        let half = (1 << VAR_BITS) / 2;
        assert!(
            structure.num_constraints <= 1 << CONSTRAINT_BITS,
            "too many constraints"
        );
        assert!(
            structure.num_io < half,
            "public part (io, 1) does not fit in the first half"
        );
        assert!(
            structure.num_witness <= half,
            "witness does not fit in the second half"
        );

        Self { a, b, c, structure }
    }

    /// io と witness から assignment $z = (io, 1, 0, \ldots \mid w, 0, \ldots)$ を組み立てる。
    pub fn assignment(&self, io: &[F], witness: &[F]) -> [F; 1 << VAR_BITS] {
        assert_eq!(io.len(), self.structure.num_io, "invalid io length");
        assert_eq!(
            witness.len(),
            self.structure.num_witness,
            "invalid witness length"
        );

        let half = (1 << VAR_BITS) / 2;
        let mut assignment = [F::zero(); 1 << VAR_BITS];
        assignment[..io.len()].copy_from_slice(io);
        assignment[io.len()] = F::one();
        assignment[half..half + witness.len()].copy_from_slice(witness);
        assignment
    }

    /// $(Az) \circ (Bz) = Cz$ を満たすかどうかを返す。
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
    use super::{Matrix, R1CS, R1CSStructure};
    use ark_bls12_381::Fr as F;

    fn f(value: u64) -> F {
        F::from(value)
    }

    /// 制約 x * y = out、witness = (x, y, out)、io なし。
    /// z = (1, 0, 0, 0 | x, y, out, 0)
    fn mul_r1cs() -> R1CS<F, 0, 3> {
        R1CS::new(
            Matrix::from_usize([[0, 0, 0, 0, 1, 0, 0, 0]]),
            Matrix::from_usize([[0, 0, 0, 0, 0, 1, 0, 0]]),
            Matrix::from_usize([[0, 0, 0, 0, 0, 0, 1, 0]]),
            R1CSStructure {
                num_constraints: 1,
                num_io: 0,
                num_witness: 3,
            },
        )
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
    fn assignment_places_public_part_and_witness_in_each_half() {
        let r1cs = mul_r1cs();

        assert_eq!(
            r1cs.assignment(&[], &[f(3), f(5), f(15)]),
            [f(1), f(0), f(0), f(0), f(3), f(5), f(15), f(0)]
        );
    }

    #[test]
    fn r1cs_accepts_satisfying_witness() {
        let r1cs = mul_r1cs();

        assert!(r1cs.is_sat(&r1cs.assignment(&[], &[f(3), f(5), f(15)])));
    }

    #[test]
    fn r1cs_rejects_unsatisfying_witness() {
        let r1cs = mul_r1cs();

        assert!(!r1cs.is_sat(&r1cs.assignment(&[], &[f(3), f(5), f(14)])));
    }
}
