//! 単項式基底の疎な多変数多項式。
//!
//! 各項を「係数 $c_i$ と変数添字の多重集合 $S_i$」で表し、
//! $$f(X_0, \ldots, X_{t-1}) = \sum_i c_i \prod_{j \in S_i} X_j$$
//! とする。多重集合なので同じ添字を重ねれば冪になる($\{0, 0\}$ は $X_0^2$)。
//! CCS の制約多項式 $\sum_i c_i \bigcirc_{j \in S_i} M_j z$ はこの $f$ に
//! $X_j = M_j z$ を(行ごとに)代入したもので、sumcheck の終点での claims
//! $v_j = \widetilde{M_j z}(r_x)$ に代入すれば $\sum_i c_i \prod_{j \in S_i} v_j$ になる。

use ark_ff::Field;

use crate::primitive::Transcript;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparseMvPoly<F: Field> {
    /// `(c_i, S_i)`。
    terms: Vec<(F, Vec<usize>)>,
}

impl<F: Field> SparseMvPoly<F> {
    pub fn new(terms: Vec<(F, Vec<usize>)>) -> Self {
        Self { terms }
    }

    /// 各項 `(c_i, S_i)`。
    pub fn terms(&self) -> &[(F, Vec<usize>)] {
        &self.terms
    }

    /// 項の数 $q$。
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// 現れる変数の添字の上限 + 1。`eval` に渡す値の最小の長さ。
    pub fn num_vars(&self) -> usize {
        self.terms
            .iter()
            .flat_map(|(_, multiset)| multiset.iter().copied())
            .max()
            .map_or(0, |j| j + 1)
    }

    /// 全次数 $\max_i |S_i|$。
    pub fn degree(&self) -> usize {
        self.terms
            .iter()
            .map(|(_, multiset)| multiset.len())
            .max()
            .unwrap_or(0)
    }

    /// $f(\vec{x}) = \sum_i c_i \prod_{j \in S_i} x_j$。
    pub fn eval(&self, point: &[F]) -> F {
        assert!(
            point.len() >= self.num_vars(),
            "point has fewer coordinates than the polynomial uses"
        );

        self.terms.iter().fold(F::zero(), |acc, (c, multiset)| {
            let monomial = multiset.iter().fold(F::one(), |prod, &j| prod * point[j]);
            acc + *c * monomial
        })
    }

    /// 項の並びごと transcript に bind する。
    pub fn append_to_transcript(&self, label: &[u8], transcript: &mut Transcript) {
        transcript.append_usize(b"mv-poly-num-terms", self.terms.len());
        for (c, multiset) in &self.terms {
            transcript.append_serializable(label, c);
            transcript.append_usize(b"mv-poly-multiset-len", multiset.len());
            for &j in multiset {
                transcript.append_usize(b"mv-poly-multiset-index", j);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::SparseMvPoly;
    use crate::primitive::Transcript;
    use ark_bls12_381::Fr as F;
    use ark_ff::One;

    /// $f = x_0^2 x_1 + x_0 + x_1 - x_2$
    fn example() -> SparseMvPoly<F> {
        SparseMvPoly::new(vec![
            (F::one(), vec![0, 0, 1]),
            (F::one(), vec![0]),
            (F::one(), vec![1]),
            (-F::one(), vec![2]),
        ])
    }

    #[test]
    fn evaluates_monomials_with_repeated_indices_as_powers() {
        let f = example();
        assert_eq!(f.degree(), 3);
        assert_eq!(f.num_vars(), 3);
        assert_eq!(f.num_terms(), 4);

        // x = (2, 3, 17): 4 * 3 + 2 + 3 - 17 = 0
        assert_eq!(f.eval(&[F::from(2), F::from(3), F::from(17)]), F::from(0));
        assert_eq!(f.eval(&[F::from(2), F::from(3), F::from(16)]), F::one());
    }

    #[test]
    fn empty_polynomial_and_empty_monomial() {
        let zero = SparseMvPoly::<F>::new(vec![]);
        assert_eq!(zero.degree(), 0);
        assert_eq!(zero.num_vars(), 0);
        assert_eq!(zero.eval(&[]), F::from(0));

        // 空の多重集合は定数項。
        let constant = SparseMvPoly::new(vec![(F::from(7), vec![])]);
        assert_eq!(constant.eval(&[F::from(3)]), F::from(7));
    }

    #[test]
    fn transcript_binding_distinguishes_term_order() {
        let f = example();
        let mut swapped_terms = f.terms().to_vec();
        swapped_terms.swap(1, 2);
        let g = SparseMvPoly::new(swapped_terms);

        let challenge = |poly: &SparseMvPoly<F>| {
            let mut transcript = Transcript::new(b"mv-poly-test");
            poly.append_to_transcript(b"poly", &mut transcript);
            transcript.challenge_field::<F>(b"probe")
        };
        assert_ne!(challenge(&f), challenge(&g));
    }
}
