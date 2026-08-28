use std::collections::BTreeMap;

use ark_ff::Field;

use super::helpers::inner_product;
use super::mle::{EqPoly, MultilinearPoly};

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
}

#[cfg(test)]
mod tests {
    use super::{DenseMatrix, Matrix, SparseMatrix};
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
        let rx = [F::from(13)];
        let ry = [F::from(17), F::from(19)];

        assert_eq!(dense.mul_vec(&vec), sparse.mul_vec(&vec));
        assert_eq!(dense.eval_mle(&rx, &ry), sparse.eval_mle(&rx, &ry));
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
}
