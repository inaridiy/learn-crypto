//! Customizable Constraint System (CCS)。
//!
//! 行列 $M_1, \ldots, M_t$、添字の多重集合 $S_1, \ldots, S_q \subseteq [t]$、
//! 定数 $c_1, \ldots, c_q$ からなり、制約は
//! $$\sum_{i=1}^{q} c_i \cdot \bigcirc_{j \in S_i} M_j z = \vec{0}$$
//! ($\bigcirc$ は Hadamard 積)。R1CS は $t = 3$, $S_1 = \{1, 2\}$, $S_2 = \{3\}$,
//! $c = (1, -1)$ の特殊ケースで、$A z \circ B z - C z = 0$ になる。
//!
//! $(S_i, c_i)$ の組は行列 $M_j z$ を変数とする多変数多項式
//! $f(X_1, \ldots, X_t) = \sum_i c_i \prod_{j \in S_i} X_j$ そのものなので、
//! [`SparseMvPoly`] として持つ。制約は行ごとに $f(M_1 z, \ldots, M_t z) = 0$。
//!
//! [`R1CS`] / [`SpartanR1CS`](crate::primitive::SpartanR1CS) と同じく、
//! dense 配置の [`CCS`] と Spartan 配置(2 冪に zero-padding)の [`SpartanCCS`] を分ける。

use std::collections::BTreeMap;
use std::marker::PhantomData;

use ark_ff::Field;

use crate::primitive::{
    DenseMultilinearPoly, Matrix, MultilinearPoly, R1CS, R1CSStructure, SparseMatrix, SparseMvPoly,
    hypercube_size,
};

/// 制約多項式 $f$ に $X_j = (M_j z)_{\mathrm{row}}$ を代入した値が、最初の `num_rows` 行で
/// すべて零か。`products[j]` は $M_j z$。
fn rows_vanish<F: Field>(poly: &SparseMvPoly<F>, products: &[Vec<F>], num_rows: usize) -> bool {
    let mut point = vec![F::zero(); products.len()];
    (0..num_rows).all(|row| {
        for (x, product) in point.iter_mut().zip(products) {
            *x = product[row];
        }
        poly.eval(&point).is_zero()
    })
}

/// 制約系 $\sum_i c_i \bigcirc_{j \in S_i} M_j z = 0$。
/// Pinocchio 流の dense な配置 `z = [1, io..., witness...]` で、padding は行わない。
pub struct CCS<F: Field, M: Matrix<F>> {
    pub matrices: Vec<M>,
    /// 制約多項式 $f(X_1, \ldots, X_t) = \sum_i c_i \prod_{j \in S_i} X_j$。
    /// 変数 $X_j$ は `matrices[j]` に対応する。
    pub poly: SparseMvPoly<F>,
    pub structure: R1CSStructure,
    _f: PhantomData<F>,
}

/// R1CS を CCS に埋め込む。$M = (A, B, C)$ で制約多項式は $X_0 X_1 - X_2$。
impl<F: Field, M: Matrix<F> + Clone> From<&R1CS<F, M>> for CCS<F, M> {
    fn from(r1cs: &R1CS<F, M>) -> Self {
        Self::new(
            vec![r1cs.a.clone(), r1cs.b.clone(), r1cs.c.clone()],
            SparseMvPoly::new(vec![(F::one(), vec![0, 1]), (-F::one(), vec![2])]),
            r1cs.structure,
        )
    }
}

impl<F: Field, M: Matrix<F>> CCS<F, M> {
    pub fn new(matrices: Vec<M>, poly: SparseMvPoly<F>, structure: R1CSStructure) -> Self {
        assert!(!matrices.is_empty(), "CCS needs at least one matrix");
        assert!(
            matrices.iter().all(|m| m.shape() == matrices[0].shape()),
            "CCS matrices must share the same shape"
        );
        assert_eq!(
            matrices[0].rows(),
            structure.num_constraints,
            "matrix rows must match the number of constraints"
        );
        assert_eq!(
            matrices[0].cols(),
            structure.dense_len(),
            "matrix columns must match the assignment length"
        );
        assert!(
            poly.num_vars() <= matrices.len(),
            "constraint polynomial refers to a matrix outside the CCS"
        );

        Self {
            matrices,
            poly,
            structure,
            _f: PhantomData,
        }
    }

    /// 行列の数 $t$。
    pub fn num_matrices(&self) -> usize {
        self.matrices.len()
    }

    /// 制約の最大次数 $d = \max_i |S_i|$。
    pub fn max_degree(&self) -> usize {
        self.poly.degree()
    }

    /// `z = [1, io..., witness...]` を作る。
    pub fn assignment(&self, io: &[F], witness: &[F]) -> Vec<F> {
        assert_eq!(io.len(), self.structure.num_io, "invalid io length");
        assert_eq!(
            witness.len(),
            self.structure.num_witness,
            "invalid witness length"
        );

        let mut z = Vec::with_capacity(self.structure.dense_len());
        z.push(F::one());
        z.extend_from_slice(io);
        z.extend_from_slice(witness);
        z
    }

    /// $f(M_1 z, \ldots, M_t z) = 0$ が全行で成り立つか。
    pub fn is_sat(&self, z: &[F]) -> bool {
        assert_eq!(
            z.len(),
            self.matrices[0].cols(),
            "invalid assignment length"
        );

        let products = self
            .matrices
            .iter()
            .map(|m| m.mul_vec(z))
            .collect::<Vec<_>>();
        rows_vanish(&self.poly, &products, self.structure.num_constraints)
    }
}

/// 片側 MLE の変数の数 $s$。$2^s \ge \max(|witness|, |io| + 1)$ を満たす最小の $s$。
fn spartan_half_vars(structure: &R1CSStructure) -> usize {
    structure
        .num_witness
        .max(structure.num_io + 1)
        .next_power_of_two()
        .ilog2() as usize
}

/// 制約側($\vec{r}_x$)の変数の数。制約行数を 2 冪に切り上げたものの $\log_2$。
fn spartan_constraint_vars(structure: &R1CSStructure) -> usize {
    structure.num_constraints.next_power_of_two().ilog2() as usize
}

/// Spartan 配置の CCS。`z = [witness..., 0..., 1, io..., 0...]` で、
/// [`SpartanR1CS`](crate::primitive::SpartanR1CS) と同じ padding 規則に従う:
/// witness 側と `(1, io)` 側をそれぞれ長さ $2^s$ に揃え、
/// $\tilde{Z}(\vec{x}, x_s) = (1 - x_s) \tilde{W}(\vec{x}) + x_s \widetilde{(1, io)}(\vec{x})$。
pub struct SpartanCCS<F: Field, M: Matrix<F>> {
    pub matrices: Vec<M>,
    /// 制約多項式。変数 $X_j$ は `matrices[j]` に対応する。
    pub poly: SparseMvPoly<F>,
    pub structure: R1CSStructure,
    _f: PhantomData<F>,
}

/// dense 配置の CCS を Spartan 配置に relax する。
/// 列を並び替え、zero-padding して 2 冪の形に揃える。
impl<F: Field, M: Matrix<F>> From<&CCS<F, M>> for SpartanCCS<F, SparseMatrix<F>> {
    fn from(ccs: &CCS<F, M>) -> Self {
        let structure = ccs.structure;
        let half_len = hypercube_size(spartan_half_vars(&structure));

        // dense 配置 [1, io..., witness...] -> Spartan 配置 [witness..., 0..., 1, io..., 0...]
        let remap = |col: usize| {
            if col == 0 {
                half_len
            } else if col <= structure.num_io {
                half_len + col
            } else {
                col - 1 - structure.num_io
            }
        };

        let matrices = ccs
            .matrices
            .iter()
            .map(|m| {
                let mut entries = BTreeMap::new();
                for (row, col, value) in m.nonzero_entries() {
                    entries.insert((row, remap(col)), value);
                }
                SparseMatrix::new(
                    entries,
                    hypercube_size(spartan_constraint_vars(&structure)),
                    2 * half_len,
                )
            })
            .collect();

        Self::new(matrices, ccs.poly.clone(), structure)
    }
}

impl<F: Field, M: Matrix<F>> SpartanCCS<F, M> {
    pub fn new(matrices: Vec<M>, poly: SparseMvPoly<F>, structure: R1CSStructure) -> Self {
        assert!(!matrices.is_empty(), "CCS needs at least one matrix");
        assert!(
            matrices.iter().all(|m| m.shape() == matrices[0].shape()),
            "CCS matrices must share the same shape"
        );
        assert_eq!(
            matrices[0].rows(),
            hypercube_size(spartan_constraint_vars(&structure)),
            "matrix rows must be the padded number of constraints"
        );
        assert_eq!(
            matrices[0].cols(),
            hypercube_size(spartan_half_vars(&structure) + 1),
            "matrix columns must match the Spartan assignment length"
        );
        assert!(
            poly.num_vars() <= matrices.len(),
            "constraint polynomial refers to a matrix outside the CCS"
        );

        Self {
            matrices,
            poly,
            structure,
            _f: PhantomData,
        }
    }

    /// 行列の数 $t$。
    pub fn num_matrices(&self) -> usize {
        self.matrices.len()
    }

    /// 制約の最大次数 $d = \max_i |S_i|$。
    pub fn max_degree(&self) -> usize {
        self.poly.degree()
    }

    /// 片側 MLE($\tilde{W}$ や $\widetilde{(1, io)}$)の変数の数 $s$。
    pub fn half_vars(&self) -> usize {
        spartan_half_vars(&self.structure)
    }

    /// $\tilde{Z}$ の変数の数 $s + 1$。
    pub fn vars(&self) -> usize {
        self.half_vars() + 1
    }

    /// 制約側($\vec{r}_x$)の変数の数。
    pub fn constraint_vars(&self) -> usize {
        spartan_constraint_vars(&self.structure)
    }

    /// 片側の長さ $2^s$。
    pub fn half_len(&self) -> usize {
        hypercube_size(self.half_vars())
    }

    /// witness 側の半分 `[witness..., 0...]` を $s$ 変数の MLE $\tilde{W}$ として返す。
    pub fn witness_mle(&self, witness: &[F]) -> DenseMultilinearPoly<F> {
        assert_eq!(
            witness.len(),
            self.structure.num_witness,
            "invalid witness length"
        );

        let mut evals = witness.to_vec();
        evals.resize(self.half_len(), F::zero());
        DenseMultilinearPoly::new(evals, self.half_vars())
    }

    /// public 側の半分 `[1, io..., 0...]` を $s$ 変数の MLE として返す。verifier が自力で計算できる。
    pub fn public_mle(&self, io: &[F]) -> DenseMultilinearPoly<F> {
        assert_eq!(io.len(), self.structure.num_io, "invalid io length");

        let mut evals = Vec::with_capacity(self.half_len());
        evals.push(F::one());
        evals.extend_from_slice(io);
        evals.resize(self.half_len(), F::zero());
        DenseMultilinearPoly::new(evals, self.half_vars())
    }

    /// `z = [witness..., 0..., 1, io..., 0...]`(長さ $2^{s+1}$)。
    pub fn assignment(&self, io: &[F], witness: &[F]) -> Vec<F> {
        let mut z = self.witness_mle(witness).into_evals();
        z.extend(self.public_mle(io).into_evals());
        z
    }

    /// $\tilde{Z}(\vec{r}, r_s) = (1 - r_s) \tilde{W}(\vec{r}) + r_s \widetilde{(1, io)}(\vec{r})$。
    pub fn assignment_eval(&self, io: &[F], point: &[F], witness_eval: F) -> F {
        assert_eq!(point.len(), self.vars(), "point dimension does not match Z");

        let (inner, selector) = (&point[..self.half_vars()], point[self.half_vars()]);
        let public_eval = self.public_mle(io).eval(inner);
        (F::one() - selector) * witness_eval + selector * public_eval
    }

    /// $f(M_1 z, \ldots, M_t z) = 0$ が最初の `num_constraints` 行で成り立つか。
    /// padding 行は $M_j z$ がすべて零なので、$f$ に定数項がない限り自明に成立する。
    pub fn is_sat(&self, z: &[F]) -> bool {
        assert_eq!(
            z.len(),
            self.matrices[0].cols(),
            "invalid assignment length"
        );

        let products = self
            .matrices
            .iter()
            .map(|m| m.mul_vec(z))
            .collect::<Vec<_>>();
        rows_vanish(&self.poly, &products, self.structure.num_constraints)
    }
}

#[cfg(test)]
mod tests {
    use super::{CCS, SpartanCCS};
    use crate::primitive::{DenseMatrix, Matrix, R1CS, R1CSStructure, SparseMvPoly};
    use ark_bls12_381::Fr as F;
    use ark_ff::One;

    /// x * y = t, t * y = out。io = [x, out], witness = [y, t] で z = [1, x, out, y, t]。
    fn square_mul_r1cs() -> R1CS<F, DenseMatrix<F, 2, 5>> {
        R1CS::new(
            DenseMatrix::from_usize([[0, 1, 0, 0, 0], [0, 0, 0, 0, 1]]),
            DenseMatrix::from_usize([[0, 0, 0, 1, 0], [0, 0, 0, 1, 0]]),
            DenseMatrix::from_usize([[0, 0, 0, 0, 1], [0, 0, 1, 0, 0]]),
            R1CSStructure::new(2, 2, 2),
        )
    }

    /// 次数 3 の CCS。io = [x, out, s], witness = [y] で z = [1, x, out, s, y]。
    ///
    /// 行 0: x * x * y + x + y - out = 0、行 1: 1 * 1 * y + 1 + y - s = 0。
    /// 制約多項式は $X_0^2 X_1 + X_0 + X_1 - X_2$。
    fn cubic_ccs() -> CCS<F, DenseMatrix<F, 2, 5>> {
        CCS::new(
            vec![
                DenseMatrix::from_usize([[0, 1, 0, 0, 0], [1, 0, 0, 0, 0]]),
                DenseMatrix::from_usize([[0, 0, 0, 0, 1], [0, 0, 0, 0, 1]]),
                DenseMatrix::from_usize([[0, 0, 1, 0, 0], [0, 0, 0, 1, 0]]),
            ],
            SparseMvPoly::new(vec![
                (F::one(), vec![0, 0, 1]),
                (F::one(), vec![0]),
                (F::one(), vec![1]),
                (-F::one(), vec![2]),
            ]),
            R1CSStructure::new(2, 3, 1),
        )
    }

    #[test]
    fn r1cs_embeds_into_ccs_with_the_same_satisfiability() {
        let r1cs = square_mul_r1cs();
        let ccs = CCS::from(&r1cs);

        assert_eq!(ccs.num_matrices(), 3);
        assert_eq!(
            ccs.poly.terms(),
            [(F::one(), vec![0, 1]), (-F::one(), vec![2])]
        );
        assert_eq!(ccs.max_degree(), 2);

        let io = [F::from(3), F::from(75)];
        let witness = [F::from(5), F::from(15)];
        assert!(r1cs.is_sat(&r1cs.assignment(&io, &witness)));
        assert!(ccs.is_sat(&ccs.assignment(&io, &witness)));

        let bad_io = [F::from(3), F::from(74)];
        assert!(!r1cs.is_sat(&r1cs.assignment(&bad_io, &witness)));
        assert!(!ccs.is_sat(&ccs.assignment(&bad_io, &witness)));
    }

    #[test]
    fn cubic_ccs_checks_high_degree_and_linear_terms() {
        let ccs = cubic_ccs();
        assert_eq!(ccs.max_degree(), 3);

        // x = 2, y = 3: out = 4 * 3 + 2 + 3 = 17, s = 3 + 1 + 3 = 7
        let io = [F::from(2), F::from(17), F::from(7)];
        let witness = [F::from(3)];
        assert!(ccs.is_sat(&ccs.assignment(&io, &witness)));

        assert!(!ccs.is_sat(&ccs.assignment(&[F::from(2), F::from(16), F::from(7)], &witness)));
        assert!(!ccs.is_sat(&ccs.assignment(&[F::from(2), F::from(17), F::from(8)], &witness)));
    }

    #[test]
    fn spartan_ccs_pads_and_remaps_like_spartan_r1cs() {
        let ccs = cubic_ccs();
        let spartan = SpartanCCS::from(&ccs);

        // half_len = max(1, 3 + 1).next_power_of_two() = 4
        assert_eq!(spartan.half_len(), 4);
        assert_eq!(spartan.half_vars(), 2);
        assert_eq!(spartan.vars(), 3);
        assert_eq!(spartan.constraint_vars(), 1);
        assert_eq!(spartan.matrices[0].shape(), (2, 8));
        assert_eq!(spartan.poly, ccs.poly);

        let io = [F::from(2), F::from(17), F::from(7)];
        let witness = [F::from(3)];
        let z = spartan.assignment(&io, &witness);
        assert_eq!(
            z,
            [3, 0, 0, 0, 1, 2, 17, 7].map(F::from),
            "witness half then public half"
        );
        assert!(spartan.is_sat(&z));
        assert!(
            !spartan.is_sat(&spartan.assignment(&[F::from(2), F::from(16), F::from(7)], &witness))
        );

        // 列の並び替えは M_j z の値を変えない。
        let dense_z = ccs.assignment(&io, &witness);
        for (dense, padded) in ccs.matrices.iter().zip(&spartan.matrices) {
            assert_eq!(dense.mul_vec(&dense_z), padded.mul_vec(&z));
        }
    }

    #[test]
    #[should_panic(expected = "constraint polynomial refers to a matrix outside the CCS")]
    fn polynomial_variables_are_validated() {
        CCS::new(
            vec![DenseMatrix::<F, 1, 2>::from_usize([[0, 1]])],
            SparseMvPoly::new(vec![(F::one(), vec![1])]),
            R1CSStructure::new(1, 1, 0),
        );
    }
}
