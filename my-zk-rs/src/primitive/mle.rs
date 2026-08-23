use std::collections::BTreeMap;

use ark_ff::Field;

use super::helpers::inner_product;

pub fn hypercube_size(vars: usize) -> usize {
    assert!(
        vars < usize::BITS as usize,
        "too many variables for usize indexing"
    );
    1usize << vars
}

/// 等式多項式 $\mathrm{eq}(\vec{r}, \vec{x}) = \prod_i (r_i x_i + (1-r_i)(1-x_i))$。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EqPoly<F: Field> {
    r: Vec<F>,
}

impl<F: Field> EqPoly<F> {
    pub fn new(r: Vec<F>) -> Self {
        Self { r }
    }

    pub fn vars(&self) -> usize {
        self.r.len()
    }

    pub fn point(&self) -> &[F] {
        &self.r
    }

    pub fn eval(&self, point: &[F]) -> F {
        assert_eq!(
            point.len(),
            self.vars(),
            "point dimension does not match polynomial"
        );

        self.r.iter().zip(point).fold(F::one(), |acc, (&r, &x)| {
            acc * (r * x + (F::one() - r) * (F::one() - x))
        })
    }

    /// `index` の bit 表現が指す Boolean な頂点での評価。
    /// $\mathrm{eq}(\vec{r}, \mathrm{bits}(i)) = \prod_j (b_j ? r_j : 1 - r_j)$
    pub fn eval_index(&self, index: usize) -> F {
        debug_assert!(index < hypercube_size(self.vars()));

        self.r.iter().enumerate().fold(F::one(), |acc, (bit, &r)| {
            acc * if (index >> bit) & 1 == 1 {
                r
            } else {
                F::one() - r
            }
        })
    }

    /// Hypercube 全点での評価表。`table()[i] == eval_index(i)` を
    /// $O(2^n)$ で一括計算する。
    pub fn table(&self) -> Vec<F> {
        let size = hypercube_size(self.vars());
        let mut result = vec![F::zero(); size];
        result[0] = F::one();

        // 変数を一つ足すごとに、既存 block を (1-r) 倍と r 倍の二つに分裂させる。
        let mut block = 1;
        for &r in &self.r {
            let one_minus_r = F::one() - r;
            let (lo, hi) = result.split_at_mut(block);

            for (lo, hi) in lo.iter_mut().zip(hi) {
                *hi = *lo * r;
                *lo *= one_minus_r;
            }

            block *= 2;
        }

        result
    }
}

pub trait MultilinearPoly<F: Field> {
    fn vars(&self) -> usize;

    /// Boolean hypercube 上の評価値を、`x_0` を最下位 bit とする順序で返す。
    fn to_evaluations(&self) -> Vec<F>;

    /// $\tilde{f}(\vec{r}) = \sum_i f(i) \cdot \mathrm{eq}(\vec{r}, \mathrm{bits}(i))$。
    fn eval(&self, point: &[F]) -> F {
        assert_eq!(
            point.len(),
            self.vars(),
            "point dimension does not match polynomial"
        );

        let weights = EqPoly::new(point.to_vec()).table();
        inner_product(&self.to_evaluations(), &weights)
    }

    /// 最下位変数 `x_0` に `r` を代入し、変数を一つ減らす。
    fn fold(&mut self, r: F);

    fn final_constant(&self) -> F {
        assert_eq!(self.vars(), 0, "polynomial is not constant yet");
        self.eval(&[])
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DenseMultilinearPoly<F: Field> {
    evals: Vec<F>,
    vars: usize,
}

impl<F: Field> DenseMultilinearPoly<F> {
    pub fn new(evals: Vec<F>, vars: usize) -> Self {
        assert_eq!(
            evals.len(),
            hypercube_size(vars),
            "evaluation table must have exactly 2^vars entries"
        );

        Self { evals, vars }
    }

    pub fn evals(&self) -> &[F] {
        &self.evals
    }

    pub fn to_sparse(&self) -> SparseMultilinearPoly<F> {
        let evals = self
            .evals
            .iter()
            .enumerate()
            .filter_map(|(index, &value)| (!value.is_zero()).then_some((index, value)))
            .collect();

        SparseMultilinearPoly::new(evals, self.vars)
    }
}

impl<F: Field> MultilinearPoly<F> for DenseMultilinearPoly<F> {
    fn vars(&self) -> usize {
        self.vars
    }

    fn to_evaluations(&self) -> Vec<F> {
        self.evals.clone()
    }

    fn fold(&mut self, r: F) {
        assert!(self.vars > 0, "cannot fold a constant polynomial");

        let one_minus_r = F::one() - r;
        let half = self.evals.len() / 2;

        for i in 0..half {
            let lo = self.evals[2 * i];
            let hi = self.evals[2 * i + 1];

            self.evals[i] = one_minus_r * lo + r * hi;
        }

        self.evals.truncate(half);
        self.vars -= 1;
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparseMultilinearPoly<F: Field> {
    evals: BTreeMap<usize, F>,
    vars: usize,
}

impl<F: Field> SparseMultilinearPoly<F> {
    pub fn new(mut evals: BTreeMap<usize, F>, vars: usize) -> Self {
        let size = hypercube_size(vars);

        assert!(
            evals.keys().all(|&index| index < size),
            "evaluation index is outside the Boolean hypercube"
        );

        evals.retain(|_, value| !value.is_zero());

        Self { evals, vars }
    }

    pub fn evals(&self) -> &BTreeMap<usize, F> {
        &self.evals
    }

    pub fn to_dense(&self) -> DenseMultilinearPoly<F> {
        DenseMultilinearPoly::new(self.to_evaluations(), self.vars)
    }
}

impl<F: Field> MultilinearPoly<F> for SparseMultilinearPoly<F> {
    fn vars(&self) -> usize {
        self.vars
    }

    fn to_evaluations(&self) -> Vec<F> {
        let mut evals = vec![F::zero(); hypercube_size(self.vars)];
        for (&index, &value) in &self.evals {
            evals[index] = value;
        }
        evals
    }

    /// 非零成分だけを走る $O(\mathrm{nnz} \cdot n)$ 版。
    fn eval(&self, point: &[F]) -> F {
        assert_eq!(
            point.len(),
            self.vars,
            "point dimension does not match polynomial"
        );

        let eq = EqPoly::new(point.to_vec());
        self.evals.iter().fold(F::zero(), |acc, (&index, &value)| {
            acc + value * eq.eval_index(index)
        })
    }

    fn fold(&mut self, r: F) {
        assert!(self.vars > 0, "cannot fold a constant polynomial");

        let one_minus_r = F::one() - r;
        let mut next = BTreeMap::new();

        for (&index, &value) in &self.evals {
            let weight = if index & 1 == 0 { one_minus_r } else { r };
            *next.entry(index >> 1).or_insert_with(F::zero) += value * weight;
        }

        // 代入や cancellation で zero になった entry は消す。
        next.retain(|_, value| !value.is_zero());

        self.evals = next;
        self.vars -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::{DenseMultilinearPoly, EqPoly, MultilinearPoly};
    use ark_bls12_381::Fr as F;

    #[test]
    fn eq_table_matches_pointwise_evaluations() {
        let eq = EqPoly::new(vec![F::from(3), F::from(5), F::from(7)]);
        let table = eq.table();

        for (index, &weight) in table.iter().enumerate() {
            assert_eq!(weight, eq.eval_index(index));
        }
    }

    #[test]
    fn dense_and_sparse_evaluations_agree() {
        let dense = DenseMultilinearPoly::new([0u64, 2, 0, 4].map(F::from).to_vec(), 2);
        let sparse = dense.to_sparse();
        let point = [F::from(11), F::from(13)];

        assert_eq!(dense.eval(&point), sparse.eval(&point));
        assert_eq!(dense.to_evaluations(), sparse.to_evaluations());
    }

    #[test]
    fn folding_evaluates_one_variable_at_a_time() {
        let mut dense = DenseMultilinearPoly::new((1..=4).map(F::from).collect(), 2);
        let mut sparse = dense.to_sparse();
        let point = [F::from(3), F::from(5)];
        let expected = dense.eval(&point);

        for &r in &point {
            dense.fold(r);
            sparse.fold(r);
        }

        assert_eq!(dense.final_constant(), expected);
        assert_eq!(sparse.final_constant(), expected);
    }
}
