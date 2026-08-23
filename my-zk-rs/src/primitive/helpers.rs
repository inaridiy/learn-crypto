use core::ops::Mul;

use ark_ff::{Field, Zero};

/// F-加群 `V` 上の内積 $\langle \vec{v}, \vec{w} \rangle = \sum_i v_i w_i$。
///
/// `V = F` なら通常の内積、`V = G`(群)なら multi-scalar multiplication になる。
pub fn inner_product<F: Field, V>(values: &[V], weights: &[F]) -> V
where
    V: Copy + Zero + Mul<F, Output = V>,
{
    assert_eq!(values.len(), weights.len(), "inner-product lengths differ");
    values
        .iter()
        .zip(weights)
        .fold(V::zero(), |sum, (&v, &w)| sum + v * w)
}

/// Split-and-fold の一段: $\vec{v} \mapsto c_L \vec{v}_L + c_R \vec{v}_R$。
///
/// `inner_product` と同じく、スカラー列にも generator 列にも使える。
pub fn fold_halves<F: Field, V>(values: &[V], left: F, right: F) -> Vec<V>
where
    V: Copy + Zero + Mul<F, Output = V>,
{
    assert_eq!(values.len() % 2, 0, "cannot halve an odd-length vector");

    let (lo, hi) = values.split_at(values.len() / 2);
    lo.iter()
        .zip(hi)
        .map(|(&l, &r)| l * left + r * right)
        .collect()
}

/// Column-major な行列の `row` 行目を列順に返す。
pub fn column_major_row<T: Copy>(
    entries: &[T],
    rows: usize,
    row: usize,
) -> impl ExactSizeIterator<Item = T> + '_ {
    assert!(rows > 0, "matrix must have at least one row");
    assert!(row < rows, "row index is outside the matrix");
    assert_eq!(
        entries.len() % rows,
        0,
        "entry count must be divisible by the row count"
    );

    entries[row..].iter().step_by(rows).copied()
}

#[cfg(test)]
mod tests {
    use super::{column_major_row, fold_halves, inner_product};
    use ark_bls12_381::Fr as F;

    #[test]
    fn extracts_a_row_from_column_major_entries() {
        let entries = [0, 1, 2, 3, 4, 5];

        assert_eq!(
            column_major_row(&entries, 2, 0).collect::<Vec<_>>(),
            [0, 2, 4]
        );
        assert_eq!(
            column_major_row(&entries, 2, 1).collect::<Vec<_>>(),
            [1, 3, 5]
        );
    }

    #[test]
    fn folds_halves_with_the_given_scales() {
        let values = [1, 2, 3, 4].map(F::from);
        let folded = fold_halves(&values, F::from(10), F::from(100));

        assert_eq!(folded, [F::from(310), F::from(420)]);
        assert_eq!(
            inner_product(&folded, &[F::from(1), F::from(1)]),
            F::from(730)
        );
    }
}
