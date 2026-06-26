use ark_ff::Field;

use super::helpers::log2_ceil;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<F: Field, const ROW_BITS: usize = 0, const COL_BITS: usize = 0>
where
    [(); 1 << ROW_BITS]:,
    [(); 1 << COL_BITS]:,
{
    entries: [[F; 1 << COL_BITS]; 1 << ROW_BITS],
}

impl<F: Field, const ROW_BITS: usize, const COL_BITS: usize> Matrix<F, ROW_BITS, COL_BITS>
where
    [(); 1 << ROW_BITS]:,
    [(); 1 << COL_BITS]:,
{
    pub fn new(entries: [[F; 1 << COL_BITS]; 1 << ROW_BITS]) -> Self {
        Self { entries }
    }

    pub fn from_usize(entries: [[usize; 1 << COL_BITS]; 1 << ROW_BITS]) -> Self {
        Self::new(entries.map(|row| row.map(|value| F::from(value as u128))))
    }

    pub fn from_unpadded<const NUM_ROWS: usize, const NUM_COLS: usize>(
        entries: [[F; NUM_COLS]; NUM_ROWS],
    ) -> Matrix<F, { log2_ceil(NUM_ROWS) }, { log2_ceil(NUM_COLS) }>
    where
        [(); 1 << log2_ceil(NUM_ROWS)]:,
        [(); 1 << log2_ceil(NUM_COLS)]:,
    {
        let mut padded = [[F::zero(); 1 << log2_ceil(NUM_COLS)]; 1 << log2_ceil(NUM_ROWS)];

        for row in 0..NUM_ROWS {
            padded[row][..NUM_COLS].copy_from_slice(&entries[row]);
        }
        Matrix::<F, { log2_ceil(NUM_ROWS) }, { log2_ceil(NUM_COLS) }>::new(padded)
    }

    pub fn from_unpadded_usize<const NUM_ROWS: usize, const NUM_COLS: usize>(
        entries: [[usize; NUM_COLS]; NUM_ROWS],
    ) -> Matrix<F, { log2_ceil(NUM_ROWS) }, { log2_ceil(NUM_COLS) }>
    where
        [(); 1 << log2_ceil(NUM_ROWS)]:,
        [(); 1 << log2_ceil(NUM_COLS)]:,
    {
        Self::from_unpadded(entries.map(|row| row.map(|value| F::from(value as u128))))
    }

    pub fn rows(&self) -> &[[F; 1 << COL_BITS]; 1 << ROW_BITS] {
        &self.entries
    }

    pub fn row(&self, index: usize) -> &[F; 1 << COL_BITS] {
        &self.entries[index]
    }

    pub fn dot_row(&self, row_index: usize, vector: &[F; 1 << COL_BITS]) -> F {
        self.entries[row_index]
            .iter()
            .zip(vector)
            .fold(F::zero(), |acc, (entry, value)| acc + (*entry * value))
    }

    pub fn mul_vector(&self, vector: &[F; 1 << COL_BITS]) -> [F; 1 << ROW_BITS] {
        std::array::from_fn(|row_index| self.dot_row(row_index, vector))
    }
}
