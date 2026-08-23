use std::collections::BTreeMap;

use ark_ff::Field;

fn hypercube_size(vars: usize) -> usize {
    debug_assert!(
        vars < usize::BITS as usize,
        "too many variables for usize indexing"
    );
    1usize << vars
}

fn set_index_to_point<F: Field>(index: usize, point: &mut [F]) {
    debug_assert!(index < hypercube_size(point.len()));

    for (i, x) in point.iter_mut().enumerate() {
        *x = if (index >> i) & 1 == 0 {
            F::zero()
        } else {
            F::one()
        };
    }
}

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

        self.r
            .iter()
            .copied()
            .zip(point.iter().copied())
            .fold(F::one(), |acc, (r, x)| {
                acc * (r * x + (F::one() - r) * (F::one() - x))
            })
    }

    pub fn table(&self) -> Vec<F> {
        let size = hypercube_size(self.vars());
        let mut point = vec![F::zero(); self.vars()];

        (0..size)
            .map(|index| {
                set_index_to_point(index, &mut point);
                self.eval(&point)
            })
            .collect()
    }
}

pub trait MultilinearPoly<F: Field> {
    fn vars(&self) -> usize;

    fn eval(&self, point: &[F]) -> F;

    fn fold(&mut self, r: F);

    fn final_constant(&self) -> F;
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
            .filter(|(_, value)| !value.is_zero())
            .map(|(index, &value)| (index, value))
            .collect();

        SparseMultilinearPoly::new(evals, self.vars)
    }
}

impl<F: Field> MultilinearPoly<F> for DenseMultilinearPoly<F> {
    fn vars(&self) -> usize {
        self.vars
    }

    fn eval(&self, point: &[F]) -> F {
        assert_eq!(
            point.len(),
            self.vars,
            "point dimension does not match polynomial"
        );

        let weights = EqPoly::new(point.to_vec()).table();

        self.evals
            .iter()
            .copied()
            .zip(weights)
            .fold(F::zero(), |acc, (eval, weight)| acc + eval * weight)
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

    fn final_constant(&self) -> F {
        assert!(self.vars() == 0);
        self.evals[0]
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
        let mut evals = vec![F::zero(); hypercube_size(self.vars)];

        for (&index, &value) in &self.evals {
            evals[index] = value;
        }

        DenseMultilinearPoly::new(evals, self.vars)
    }
}

impl<F: Field> MultilinearPoly<F> for SparseMultilinearPoly<F> {
    fn vars(&self) -> usize {
        self.vars
    }

    fn eval(&self, point: &[F]) -> F {
        assert_eq!(
            point.len(),
            self.vars,
            "point dimension does not match polynomial"
        );

        let eq = EqPoly::new(point.to_vec());
        let mut boolean_point = vec![F::zero(); self.vars];

        self.evals.iter().fold(F::zero(), |acc, (&index, &eval)| {
            set_index_to_point(index, &mut boolean_point);
            acc + eval * eq.eval(&boolean_point)
        })
    }

    fn fold(&mut self, r: F) {
        assert!(self.vars > 0, "cannot fold a constant polynomial");

        let one_minus_r = F::one() - r;
        let mut next = BTreeMap::new();

        for (&index, &eval) in &self.evals {
            let weight = if index & 1 == 0 { one_minus_r } else { r };

            let contribution = eval * weight;

            if contribution.is_zero() {
                continue;
            }

            *next.entry(index >> 1).or_insert_with(F::zero) += contribution;
        }

        // cancellation により zero になった entry も消す。
        next.retain(|_, value| !value.is_zero());

        self.evals = next;
        self.vars -= 1;
    }

    fn final_constant(&self) -> F {
        assert_eq!(self.vars(), 0, "polynomial is not constant yet");
        self.eval(&[])
    }
}
