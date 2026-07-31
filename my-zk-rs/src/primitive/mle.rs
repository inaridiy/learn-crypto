use ark_ff::Field;
use std::ops::Index;

use crate::primitive::MvPolynomial;

use super::helpers::log2_ceil;
use super::matrix::Matrix;

/// `N` 個の真偽値からなるブール点。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BoolPoint<const N: usize> {
    coordinates: [bool; N],
}

impl<const N: usize> BoolPoint<N> {
    /// 真偽値配列からブール点を作る。
    #[inline]
    pub fn new(coordinates: [bool; N]) -> Self {
        Self { coordinates }
    }

    /// 保持している座標を返す。
    #[inline]
    pub fn coordinates(&self) -> &[bool; N] {
        &self.coordinates
    }

    /// テーブル添字に変換する。
    ///
    /// `x0` を最下位ビット、`x{N-1}` を最上位ビットとして扱う。
    #[inline]
    pub fn to_index(&self) -> usize {
        self.coordinates
            .iter()
            .enumerate()
            .fold(0, |index, (i, bit)| index | ((*bit as usize) << i))
    }

    #[inline]
    pub fn to_field_point<F: Field>(&self) -> [F; N] {
        std::array::from_fn(|i| {
            if self.coordinates[i] {
                F::one()
            } else {
                F::zero()
            }
        })
    }

    /// このブール点に対応する equality polynomial を返す。
    ///
    /// 戻り値 `eq_w(x)` は、ブール点上で `x == w` のときだけ `1`、それ以外で `0` になる。
    #[inline]
    pub fn teq<F: Field>(&self) -> MvPolynomial<F, N> {
        let mut eq = MvPolynomial::one();

        for i in 0..N {
            let variable = MvPolynomial::variable(i);
            let factor = if self.coordinates[i] {
                variable
            } else {
                MvPolynomial::one() - variable
            };
            eq *= factor;
        }

        eq
    }

    #[inline]
    pub fn eq_eval<F: Field>(&self, point: &[F; N]) -> F {
        self.coordinates
            .iter()
            .enumerate()
            .fold(F::one(), |acc, (i, bit)| {
                if *bit {
                    acc * point[i]
                } else {
                    acc * (F::one() - point[i])
                }
            })
    }
}

impl<const N: usize> From<[bool; N]> for BoolPoint<N> {
    #[inline]
    fn from(value: [bool; N]) -> Self {
        Self::new(value)
    }
}

impl<const N: usize> From<usize> for BoolPoint<N> {
    #[inline]
    fn from(value: usize) -> Self {
        let coordinates = std::array::from_fn(|i| {
            if i < usize::BITS as usize {
                ((value >> i) & 1) == 1
            } else {
                false
            }
        });
        Self { coordinates }
    }
}

impl<const N: usize> Index<usize> for BoolPoint<N> {
    type Output = bool;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.coordinates[index]
    }
}

/// ブール超立方体 `{0, 1}^N`。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BoolHyperCube<const N: usize>
where
    [(); 1 << N]:;

impl<const N: usize> BoolHyperCube<N>
where
    [(); 1 << N]:,
{
    /// ブール超立方体上の点をテーブル添字順に列挙する。
    #[inline]
    pub fn iter() -> impl ExactSizeIterator<Item = BoolPoint<N>> {
        (0..Self::len()).map(BoolPoint::<N>::from)
    }

    /// ブール超立方体に含まれる点の個数。
    #[inline]
    pub const fn len() -> usize {
        1 << N
    }

    /// ブール超立方体は常に少なくとも 1 点を持つ。
    #[inline]
    pub const fn is_empty() -> bool {
        false
    }

    /// 各ブール点 `w` について equality polynomial `eq_w(point)` をテーブル順に並べる。
    #[inline]
    pub fn eq_evaluations<F: Field>(point: &[F; N]) -> [F; 1 << N] {
        std::array::from_fn(|i| BoolPoint::<N>::from(i).eq_eval(point))
    }
}

/// Multilinear extension の dense evaluation table。
///
/// `evaluations[point.to_index()]` がブール点 `point` での関数値を表す。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DenseFnOverBoolHyperCube<F: Field, const N: usize>
where
    [(); 1 << N]:,
{
    evaluations: [F; 1 << N],
}

impl<F: Field, const N: usize> DenseFnOverBoolHyperCube<F, N>
where
    [(); 1 << N]:,
{
    /// ブール超立方体上の評価値から MLE を作る。
    #[inline]
    pub fn new(evaluations: [F; 1 << N]) -> Self {
        Self { evaluations }
    }

    /// ブール超立方体上の評価値を返す。
    #[inline]
    pub fn evaluations(&self) -> &[F; 1 << N] {
        &self.evaluations
    }

    /// ブール点で元の関数を評価する。
    #[inline]
    pub fn eval(&self, point: BoolPoint<N>) -> F {
        self.evaluations[point.to_index()]
    }

    /// 定義通りに `sum_w f(w) * eq_w(x)` として MLE を構成する。

    #[inline]
    pub fn mle(&self) -> MvPolynomial<F, N> {
        let mut polynomial = MvPolynomial::zero();

        for point in BoolHyperCube::<N>::iter() {
            let evaluation = self.eval(point);
            if evaluation.is_zero() {
                continue;
            }

            polynomial += point.teq::<F>().scale(evaluation);
        }

        polynomial
    }
}

impl<F: Field, const N: usize> From<[F; 1 << N]> for DenseFnOverBoolHyperCube<F, N>
where
    [(); 1 << N]:,
{
    #[inline]
    fn from(evaluations: [F; 1 << N]) -> Self {
        Self::new(evaluations)
    }
}

/// `{0, 1}^N` 上の評価値列から multilinear extension を作る。
#[inline]
pub fn mle_from_hypercube_evaluations<F: Field, const N: usize>(
    evaluations: [F; 1 << N],
) -> MvPolynomial<F, N>
where
    [(); 1 << N]:,
{
    DenseFnOverBoolHyperCube::<F, N>::from(evaluations).mle()
}

/// 評価値列から multilinear extension を作る。
///
/// 入力長 `M` に対して `M <= 2^n` を満たす最小の `n` を const 計算し、
/// 足りない評価値は `0` で埋めた上で `{0, 1}^n` 上の dense table として扱う。
#[inline]
pub fn mle_from_evaluations<F: Field, const M: usize>(
    evaluations: [F; M],
) -> MvPolynomial<F, { log2_ceil(M) }>
where
    [(); 1 << log2_ceil(M)]:,
{
    let padded_evaluations: [F; 1 << log2_ceil(M)] =
        std::array::from_fn(|i| if i < M { evaluations[i] } else { F::zero() });

    DenseFnOverBoolHyperCube::<F, { log2_ceil(M) }>::from(padded_evaluations).mle()
}

/// Matrix を `(row, col)` の MLE にする。
#[inline]
pub fn mle_from_matrix<F, const ROW_BITS: usize, const COL_BITS: usize>(
    matrix: &Matrix<F, ROW_BITS, COL_BITS>,
) -> MvPolynomial<F, { ROW_BITS + COL_BITS }>
where
    F: Field,
    [(); 1 << ROW_BITS]:,
    [(); 1 << COL_BITS]:,
    [(); 1 << (ROW_BITS + COL_BITS)]:,
{
    let evaluations: [F; 1 << (ROW_BITS + COL_BITS)] = std::array::from_fn(|i| {
        let row = i & ((1 << ROW_BITS) - 1);
        let col = i >> ROW_BITS;
        matrix.row(row)[col]
    });

    DenseFnOverBoolHyperCube::<F, { ROW_BITS + COL_BITS }>::new(evaluations).mle()
}

#[cfg(test)]
mod tests {
    use super::{BoolHyperCube, BoolPoint, DenseFnOverBoolHyperCube, mle_from_evaluations};
    use crate::primitive::Matrix;
    use ark_bls12_381::Fr as F;
    use ark_ff::{AdditiveGroup, Field};

    fn f(x: u64) -> F {
        F::from(x)
    }

    #[test]
    fn bool_point_converts_to_index_and_field_point() {
        let point = BoolPoint::<3>::new([true, false, true]);

        assert_eq!(point.to_index(), 5);
        assert_eq!(point.to_field_point::<F>(), [F::ONE, F::ZERO, F::ONE]);
        assert!(point[0]);
        assert!(!point[1]);
    }

    #[test]
    fn bool_hypercube_iterates_in_table_order() {
        let points = BoolHyperCube::<2>::iter().collect::<Vec<_>>();

        assert_eq!(
            points,
            vec![
                BoolPoint::from([false, false]),
                BoolPoint::from([true, false]),
                BoolPoint::from([false, true]),
                BoolPoint::from([true, true]),
            ]
        );
        assert_eq!(BoolHyperCube::<2>::iter().len(), 4);
        assert!(!BoolHyperCube::<2>::is_empty());
    }

    #[test]
    fn dense_mle_evaluates_on_boolean_points() {
        let mle = DenseFnOverBoolHyperCube::<F, 2>::from([f(3), f(5), f(7), f(11)]);

        assert_eq!(mle.eval([false, false].into()), f(3));
        assert_eq!(mle.eval([true, false].into()), f(5));
        assert_eq!(mle.eval([false, true].into()), f(7));
        assert_eq!(mle.eval([true, true].into()), f(11));
    }

    #[test]
    fn teq_is_one_only_at_matching_boolean_point() {
        let point = BoolPoint::<2>::from([true, false]);
        let teq = point.teq::<F>();

        assert_eq!(
            teq.eval(&BoolPoint::from([false, false]).to_field_point()),
            F::ZERO
        );
        assert_eq!(
            teq.eval(&BoolPoint::from([true, false]).to_field_point()),
            F::ONE
        );
        assert_eq!(
            teq.eval(&BoolPoint::from([false, true]).to_field_point()),
            F::ZERO
        );
        assert_eq!(
            teq.eval(&BoolPoint::from([true, true]).to_field_point()),
            F::ZERO
        );
    }

    #[test]
    fn eq_eval_matches_teq_eval() {
        let point = [f(2), f(3), f(5)];

        for bool_point in BoolHyperCube::<3>::iter() {
            assert_eq!(
                bool_point.eq_eval(&point),
                bool_point.teq::<F>().eval(&point)
            );
        }
    }

    #[test]
    fn bool_hypercube_eq_evaluations_are_in_table_order() {
        let point = [f(2), f(3)];
        let evaluations = BoolHyperCube::<2>::eq_evaluations(&point);

        assert_eq!(
            evaluations,
            [
                (F::ONE - point[0]) * (F::ONE - point[1]),
                point[0] * (F::ONE - point[1]),
                (F::ONE - point[0]) * point[1],
                point[0] * point[1],
            ]
        );
    }

    #[test]
    fn mle_matches_dense_evaluations_on_boolean_points() {
        let table = DenseFnOverBoolHyperCube::<F, 2>::from([f(3), f(5), f(7), f(11)]);
        let mle = table.mle();

        for point in BoolHyperCube::<2>::iter() {
            assert_eq!(mle.eval(&point.to_field_point()), table.eval(point));
        }
        assert_eq!(mle.degrees(), [1, 1]);
    }

    #[test]
    fn mle_from_evaluations_pads_missing_boolean_points_with_zero() {
        let mle = mle_from_evaluations([f(3), f(5), f(7)]);

        assert_eq!(mle.eval(&[f(0), f(0)]), f(3));
        assert_eq!(mle.eval(&[f(1), f(0)]), f(5));
        assert_eq!(mle.eval(&[f(0), f(1)]), f(7));
        assert_eq!(mle.eval(&[f(1), f(1)]), f(0));
    }

    #[test]
    fn mle_from_evaluations_accepts_an_empty_table_as_zero_constant() {
        let mle = mle_from_evaluations::<F, 0>([]);

        assert_eq!(mle.eval(&[]), f(0));
        assert!(mle.is_zero());
    }
}
