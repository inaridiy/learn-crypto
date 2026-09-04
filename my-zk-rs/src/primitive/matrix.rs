use std::collections::BTreeMap;

use ark_ff::Field;

use super::helpers::inner_product;
use super::mle::EqPoly;

pub trait Matrix<F: Field> {
    fn rows(&self) -> usize;
    fn cols(&self) -> usize;
    fn shape(&self) -> (usize, usize) {
        (self.rows(), self.cols())
    }

    fn get(&self, row: usize, col: usize) -> F;

    fn row(&self, row: usize) -> impl Iterator<Item = F>;

    /// $M \vec{v}$。
    fn mul_vec(&self, vec: &[F]) -> Vec<F> {
        assert_eq!(
            vec.len(),
            self.cols(),
            "vector length does not match matrix columns"
        );

        (0..self.rows())
            .map(|row| {
                self.row(row)
                    .zip(vec)
                    .fold(F::zero(), |sum, (m, &v)| sum + m * v)
            })
            .collect()
    }

    /// $\vec{v}^\top M$。行方向を潰すので `mul_vec` の転置版。
    fn vec_mul(&self, vec: &[F]) -> Vec<F> {
        assert_eq!(
            vec.len(),
            self.rows(),
            "vector length does not match matrix rows"
        );

        (0..self.cols())
            .map(|col| {
                (0..self.rows()).fold(F::zero(), |sum, row| sum + vec[row] * self.get(row, col))
            })
            .collect()
    }

    /// Multilinear extension の評価
    /// $\tilde{M}(r_x, r_y) = \mathrm{eq}(r_x)^\top M \, \mathrm{eq}(r_y)$。
    fn eval_mle(&self, rx: &[F], ry: &[F]) -> F {
        let row_weights = EqPoly::new(rx.to_vec()).to_evals();
        let col_weights = EqPoly::new(ry.to_vec()).to_evals();

        inner_product(&self.mul_vec(&col_weights), &row_weights)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DenseMatrix<F: Field, const ROWS: usize, const COLS: usize> {
    entries: [[F; COLS]; ROWS],
}

impl<F: Field, const ROWS: usize, const COLS: usize> DenseMatrix<F, ROWS, COLS> {
    pub fn new(entries: [[F; COLS]; ROWS]) -> Self {
        Self { entries }
    }

    pub fn from_usize(entries: [[usize; COLS]; ROWS]) -> Self {
        Self::new(entries.map(|row| row.map(|value| F::from(value as u128))))
    }
}

impl<F: Field, const ROWS: usize, const COLS: usize> Matrix<F> for DenseMatrix<F, ROWS, COLS> {
    fn rows(&self) -> usize {
        ROWS
    }

    fn cols(&self) -> usize {
        COLS
    }

    fn get(&self, row: usize, col: usize) -> F {
        self.entries[row][col]
    }

    fn row(&self, row: usize) -> impl Iterator<Item = F> {
        self.entries[row].iter().copied()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparseMatrix<F: Field> {
    entries: BTreeMap<(usize, usize), F>,
    rows: usize,
    cols: usize,
}

impl<F: Field> SparseMatrix<F> {
    pub fn new(mut entries: BTreeMap<(usize, usize), F>, rows: usize, cols: usize) -> Self {
        assert!(
            entries.keys().all(|&(row, col)| row < rows && col < cols),
            "entry index is outside the matrix"
        );

        entries.retain(|_, value| !value.is_zero());

        Self {
            entries,
            rows,
            cols,
        }
    }
}

impl<F: Field> Matrix<F> for SparseMatrix<F> {
    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn get(&self, row: usize, col: usize) -> F {
        self.entries.get(&(row, col)).copied().unwrap_or_default()
    }

    fn row(&self, row: usize) -> impl Iterator<Item = F> {
        (0..self.cols).map(move |col| self.get(row, col))
    }

    /// 非零成分だけを走る $O(\mathrm{nnz})$ 版。
    fn mul_vec(&self, vec: &[F]) -> Vec<F> {
        assert_eq!(
            vec.len(),
            self.cols,
            "vector length does not match matrix columns"
        );

        let mut out = vec![F::zero(); self.rows];
        for (&(row, col), &value) in &self.entries {
            out[row] += value * vec[col];
        }
        out
    }

    /// 非零成分だけを走る $O(\mathrm{nnz})$ 版。
    fn vec_mul(&self, vec: &[F]) -> Vec<F> {
        assert_eq!(
            vec.len(),
            self.rows,
            "vector length does not match matrix rows"
        );

        let mut out = vec![F::zero(); self.cols];
        for (&(row, col), &value) in &self.entries {
            out[col] += vec[row] * value;
        }
        out
    }
}

/// 一次元の評価表を column-major な `rows` 行の行列と見なす view。
/// `get(row, col) = entries[row + rows * col]` で、下位 bit が行、上位 bit が列を選ぶ。
/// Hyrax が MLE の評価表を行列 `T` として扱うときの配置。
pub struct ColumnMajorMatrix<F: Field> {
    entries: Vec<F>,
    rows: usize,
}

impl<F: Field> ColumnMajorMatrix<F> {
    pub fn new(entries: Vec<F>, rows: usize) -> Self {
        assert!(rows > 0, "matrix must have at least one row");
        assert_eq!(
            entries.len() % rows,
            0,
            "entry count must be divisible by the row count"
        );

        Self { entries, rows }
    }
}

impl<F: Field> Matrix<F> for ColumnMajorMatrix<F> {
    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.entries.len() / self.rows
    }

    fn get(&self, row: usize, col: usize) -> F {
        self.entries[row + self.rows * col]
    }

    fn row(&self, row: usize) -> impl Iterator<Item = F> {
        assert!(row < self.rows, "row index is outside the matrix");
        self.entries[row..].iter().step_by(self.rows).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::{ColumnMajorMatrix, DenseMatrix, Matrix, SparseMatrix};
    use crate::primitive::{DenseMultilinearPoly, MultilinearPoly};
    use ark_bls12_381::Fr as F;
    use std::collections::BTreeMap;

    fn sparse_from_dense<const R: usize, const C: usize>(
        dense: &DenseMatrix<F, R, C>,
    ) -> SparseMatrix<F> {
        let mut entries = BTreeMap::new();
        for row in 0..R {
            for col in 0..C {
                entries.insert((row, col), dense.get(row, col));
            }
        }
        SparseMatrix::new(entries, R, C)
    }

    #[test]
    fn dense_and_sparse_matrices_agree_on_mul_vec_and_mle() {
        let dense = DenseMatrix::<F, 2, 4>::from_usize([[1, 2, 0, 4], [0, 6, 7, 0]]);
        let sparse = sparse_from_dense(&dense);
        let vec = [3, 5, 7, 11].map(F::from);
        let row_vec = [23, 29].map(F::from);
        let rx = [F::from(13)];
        let ry = [F::from(17), F::from(19)];

        assert_eq!(dense.mul_vec(&vec), sparse.mul_vec(&vec));
        assert_eq!(dense.vec_mul(&row_vec), sparse.vec_mul(&row_vec));
        assert_eq!(dense.eval_mle(&rx, &ry), sparse.eval_mle(&rx, &ry));
    }

    #[test]
    fn vec_mul_is_the_transpose_of_mul_vec() {
        let dense = DenseMatrix::<F, 2, 3>::from_usize([[1, 2, 3], [4, 5, 6]]);
        let row_vec = [7, 11].map(F::from);

        // v^T M = (7*1 + 11*4, 7*2 + 11*5, 7*3 + 11*6)
        assert_eq!(dense.vec_mul(&row_vec), [51, 69, 87].map(F::from).to_vec());
    }

    #[test]
    fn mle_extends_the_matrix_entries() {
        let dense = DenseMatrix::<F, 2, 2>::from_usize([[1, 2], [3, 4]]);

        // Boolean な点では行列の成分そのもの。
        for row in 0..2usize {
            for col in 0..2usize {
                assert_eq!(
                    dense.eval_mle(&[F::from(row as u64)], &[F::from(col as u64)]),
                    dense.get(row, col)
                );
            }
        }
    }

    #[test]
    fn column_major_matrix_reads_rows_with_a_stride() {
        let m = ColumnMajorMatrix::new([0u64, 1, 2, 3, 4, 5].map(F::from).to_vec(), 2);

        assert_eq!(m.shape(), (2, 3));
        assert_eq!(m.row(0).collect::<Vec<_>>(), [0u64, 2, 4].map(F::from));
        assert_eq!(m.row(1).collect::<Vec<_>>(), [1u64, 3, 5].map(F::from));
        assert_eq!(m.get(1, 2), F::from(5));
    }

    #[test]
    fn column_major_mle_evaluation_splits_into_row_and_column_points() {
        // 3 変数の MLE を 2 行 4 列と見なすと、下位 1 変数が行、上位 2 変数が列。
        let poly = DenseMultilinearPoly::new((1..=8).map(F::from).collect(), 3);
        let t = ColumnMajorMatrix::new(poly.evals().to_vec(), 2);
        let point = [F::from(2), F::from(3), F::from(5)];

        assert_eq!(t.eval_mle(&point[..1], &point[1..]), poly.eval(&point));
    }
}
