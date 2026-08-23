use ark_ff::Field;

pub trait Matrix<T> {
    fn rows(&self) -> usize;
    fn cols(&self) -> usize;
    fn get(&self, row: usize, col: usize) -> &T;

    fn row<'a>(&'a self, row: usize) -> impl Iterator<Item = &'a T> + 'a
    where
        T: 'a;
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

    fn get(&self, row: usize, col: usize) -> &F {
        &self.entries[row][col]
    }

    fn row<'a>(&self, row: usize) -> impl Iterator<Item = &F> + '_ {
        self.entries[row].iter()
    }
}
