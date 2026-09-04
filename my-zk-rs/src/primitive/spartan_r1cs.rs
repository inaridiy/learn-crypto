use std::collections::BTreeMap;
use std::marker::PhantomData;

use ark_ff::Field;

use crate::primitive::{
    hypercube_size, DenseMultilinearPoly, Matrix, MultilinearPoly, R1CSStructure, SparseMatrix,
    R1CS,
};

/// Spartan 流の制約系。`z = [witness..., 0..., 1, io..., 0...]` で、
/// witness 側と `(1, io)` 側をそれぞれ長さ $2^s$ に zero-padding し全長を $2^{s+1}$ にする。
/// この repo の MLE は `x_0` が最下位 bit なので、最上位変数 $x_s$ が半分を選び
/// $\tilde{Z}(\vec{x}, x_s) = (1 - x_s) \tilde{W}(\vec{x}) + x_s \widetilde{(1, io)}(\vec{x})$。
pub struct SpartanR1CS<F: Field, M: Matrix<F>> {
    pub a: M,
    pub b: M,
    pub c: M,
    pub structure: R1CSStructure,
    _f: PhantomData<F>,
}

/// 純粋な R1CS を Spartan 配置に relax する。
/// dense 配置 `[1, io..., witness...]` の列を並び替え、zero-padding して 2 冪の形に揃える。
impl<F: Field, M: Matrix<F>> From<&R1CS<F, M>> for SpartanR1CS<F, SparseMatrix<F>> {
    fn from(r1cs: &R1CS<F, M>) -> Self {
        let structure = r1cs.structure;
        let half_len = hypercube_size(Self::half_vars_of(&structure));

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

        let convert = |m: &M| {
            let mut entries = BTreeMap::new();
            for (row, col, value) in m.nonzero_entries() {
                entries.insert((row, remap(col)), value);
            }
            SparseMatrix::new(
                entries,
                hypercube_size(Self::constraint_vars_of(&structure)),
                2 * half_len,
            )
        };

        Self::new(
            convert(&r1cs.a),
            convert(&r1cs.b),
            convert(&r1cs.c),
            structure,
        )
    }
}

impl<F: Field, M: Matrix<F>> SpartanR1CS<F, M> {
    pub fn new(a: M, b: M, c: M, structure: R1CSStructure) -> Self {
        assert!(
            a.shape() == b.shape() && a.shape() == c.shape(),
            "R1CS matrices must share the same shape"
        );
        assert_eq!(
            a.rows(),
            hypercube_size(Self::constraint_vars_of(&structure)),
            "matrix rows must be the padded number of constraints"
        );
        assert_eq!(
            a.cols(),
            hypercube_size(Self::half_vars_of(&structure) + 1),
            "matrix columns must match the Spartan assignment length"
        );

        Self {
            a,
            b,
            c,
            structure,
            _f: PhantomData,
        }
    }

    /// 片側 MLE の変数の数 $s$。$2^s \ge \max(|witness|, |io| + 1)$ を満たす最小の $s$。
    /// witness 側と `(1, io)` 側を同じ長さに揃えることで、
    /// $\tilde{Z}$ を二つの MLE の線形結合に分解できる。
    fn half_vars_of(structure: &R1CSStructure) -> usize {
        structure
            .num_witness
            .max(structure.num_io + 1)
            .next_power_of_two()
            .ilog2() as usize
    }

    /// 制約側($\vec{r}_x$)の変数の数。制約行数を 2 冪に切り上げたものの $\log_2$。
    fn constraint_vars_of(structure: &R1CSStructure) -> usize {
        structure.num_constraints.next_power_of_two().ilog2() as usize
    }

    /// 片側 MLE($\tilde{W}$ や $\widetilde{(1, io)}$)の変数の数 $s$。
    pub fn half_vars(&self) -> usize {
        Self::half_vars_of(&self.structure)
    }

    /// $\tilde{Z}$ の変数の数 $s + 1$。
    pub fn vars(&self) -> usize {
        self.half_vars() + 1
    }

    /// 制約側($\vec{r}_x$)の変数の数。
    pub fn constraint_vars(&self) -> usize {
        Self::constraint_vars_of(&self.structure)
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
    /// $\tilde{Z}$ として使うときは [`DenseMultilinearPoly::from_evals`] に渡す。
    pub fn assignment(&self, io: &[F], witness: &[F]) -> Vec<F> {
        let mut z = self.witness_mle(witness).into_evals();
        z.extend(self.public_mle(io).into_evals());
        z
    }

    /// $\tilde{Z}(\vec{r}, r_s) = (1 - r_s) \tilde{W}(\vec{r}) + r_s \widetilde{(1, io)}(\vec{r})$。
    /// verifier は $\tilde{W}(\vec{r})$ を prover から受け取り、public 側は自分で評価する。
    pub fn assignment_eval(&self, io: &[F], point: &[F], witness_eval: F) -> F {
        assert_eq!(point.len(), self.vars(), "point dimension does not match Z");

        let (inner, selector) = (&point[..self.half_vars()], point[self.half_vars()]);
        let public_eval = self.public_mle(io).eval(inner);
        (F::one() - selector) * witness_eval + selector * public_eval
    }

    /// $A z \circ B z = C z$ が最初の `num_constraints` 行で成り立つか。
    /// padding 行はすべて零なので自明に成立する。
    pub fn is_sat(&self, z: &[F]) -> bool {
        assert_eq!(z.len(), self.a.cols(), "invalid assignment length");

        let az = self.a.mul_vec(z);
        let bz = self.b.mul_vec(z);
        let cz = self.c.mul_vec(z);

        az.iter()
            .zip(&bz)
            .zip(&cz)
            .take(self.structure.num_constraints)
            .all(|((&a, &b), &c)| a * b == c)
    }
}

#[cfg(test)]
mod tests {
    use super::SpartanR1CS;
    use crate::primitive::{
        DenseMatrix, DenseMultilinearPoly, Matrix, MultilinearPoly, R1CSStructure, R1CS,
    };
    use ark_bls12_381::Fr as F;
    use ark_ff::{One, Zero};

    /// x * y = out。io = [x, out], witness = [y] で z = [1, x, out, y]。
    fn mul_r1cs() -> R1CS<F, DenseMatrix<F, 1, 4>> {
        R1CS::new(
            DenseMatrix::from_usize([[0, 1, 0, 0]]),
            DenseMatrix::from_usize([[0, 0, 0, 1]]),
            DenseMatrix::from_usize([[0, 0, 1, 0]]),
            R1CSStructure::new(1, 2, 1),
        )
    }

    #[test]
    fn spartan_assignment_places_witness_and_public_halves() {
        let spartan = SpartanR1CS::from(&mul_r1cs());

        // half_len = max(1, 2 + 1).next_power_of_two() = 4
        assert_eq!(spartan.half_len(), 4);
        assert_eq!(spartan.a.shape(), (1, 8));
        assert_eq!(spartan.half_vars(), 2);
        assert_eq!(spartan.vars(), 3);
        assert_eq!(spartan.constraint_vars(), 0);

        let z = spartan.assignment(&[F::from(3), F::from(15)], &[F::from(5)]);
        assert_eq!(
            z,
            [
                // witness 側: [y, 0, 0, 0]
                F::from(5),
                F::zero(),
                F::zero(),
                F::zero(),
                // public 側: [1, x, out, 0]
                F::one(),
                F::from(3),
                F::from(15),
                F::zero(),
            ]
        );
    }

    #[test]
    fn spartan_satisfaction_agrees_with_dense() {
        let r1cs = mul_r1cs();
        let spartan = SpartanR1CS::from(&r1cs);

        let io = [F::from(3), F::from(15)];
        let witness = [F::from(5)];
        assert!(r1cs.is_sat(&r1cs.assignment(&io, &witness)));
        assert!(spartan.is_sat(&spartan.assignment(&io, &witness)));

        let bad_io = [F::from(3), F::from(14)];
        assert!(!r1cs.is_sat(&r1cs.assignment(&bad_io, &witness)));
        assert!(!spartan.is_sat(&spartan.assignment(&bad_io, &witness)));
    }

    #[test]
    fn assignment_eval_splits_into_witness_and_public_halves() {
        let spartan = SpartanR1CS::from(&mul_r1cs());
        let io = [F::from(3), F::from(15)];
        let witness = [F::from(5)];

        let point = [F::from(7), F::from(11), F::from(13)];
        let expected =
            DenseMultilinearPoly::from_evals(spartan.assignment(&io, &witness)).eval(&point);
        let witness_eval = spartan.witness_mle(&witness).eval(&point[..2]);

        assert_eq!(spartan.assignment_eval(&io, &point, witness_eval), expected);
    }
}
