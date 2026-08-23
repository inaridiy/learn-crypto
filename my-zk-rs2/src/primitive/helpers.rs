use ark_ec::CurveGroup;
use ark_ff::Field;

use super::{Monomial, MvPolynomial};

/// Create the additive identity with the same type as `self`.
pub trait ZeroLike {
    fn zero_like(&self) -> Self;
}

/// Create the multiplicative identity with the same type as `self`.
pub trait OneLike {
    fn one_like(&self) -> Self;
}

/// Create a constant value with the same ambient type as `self`.
pub trait ConstantLike<Coeff> {
    fn constant_like(&self, coeff: Coeff) -> Self;
}

/// Create a variable with the same ambient type as `self`.
pub trait VariableLike {
    fn variable_like(&self, index: usize) -> Self;
}

impl<const N: usize> OneLike for Monomial<N> {
    #[inline]
    fn one_like(&self) -> Self {
        Self::one()
    }
}

impl<const N: usize> VariableLike for Monomial<N> {
    #[inline]
    fn variable_like(&self, index: usize) -> Self {
        Self::variable(index)
    }
}

impl<F: Field, const N: usize> ZeroLike for MvPolynomial<F, N> {
    #[inline]
    fn zero_like(&self) -> Self {
        Self::zero()
    }
}

impl<F: Field, const N: usize> OneLike for MvPolynomial<F, N> {
    #[inline]
    fn one_like(&self) -> Self {
        Self::one()
    }
}

impl<F: Field, const N: usize> ConstantLike<F> for MvPolynomial<F, N> {
    #[inline]
    fn constant_like(&self, coeff: F) -> Self {
        Self::constant(coeff)
    }
}

impl<F: Field, const N: usize> VariableLike for MvPolynomial<F, N> {
    #[inline]
    fn variable_like(&self, index: usize) -> Self {
        Self::variable(index)
    }
}

/// `value_count <= 2^n` を満たす最小の `n` を返す。
#[inline]
pub const fn log2_ceil(value_count: usize) -> usize {
    if value_count <= 1 {
        return 0;
    }

    usize::BITS as usize - (value_count - 1).leading_zeros() as usize
}

#[inline]
pub fn inner_product<F: Field>(lhs: &[F], rhs: &[F]) -> F {
    assert_eq!(lhs.len(), rhs.len(), "inner-product vector lengths differ");
    lhs.iter()
        .zip(rhs)
        .fold(F::zero(), |acc, (lhs, rhs)| acc + *lhs * rhs)
}

#[inline]
pub fn msm_with_bases<G: CurveGroup>(values: &[G::ScalarField], bases: &[G::Affine]) -> G {
    G::msm_unchecked(bases, values)
}

#[inline]
pub fn msm<G: CurveGroup>(values: &[G::ScalarField], generators: &[G]) -> G {
    let bases = G::batch_convert_to_mul_base(generators);
    msm_with_bases(values, &bases)
}

#[cfg(test)]
mod tests {
    use super::{ConstantLike, OneLike, VariableLike, ZeroLike, inner_product, log2_ceil};
    use crate::primitive::{Monomial, MvPolynomial};
    use ark_bls12_381::Fr as F;

    type Poly<const N: usize> = MvPolynomial<F, N>;

    #[test]
    fn monomial_helpers_match_receiver_dimension() {
        let x0 = Monomial::<3>::variable(0);

        assert_eq!(x0.one_like(), Monomial::<3>::one());
        assert_eq!(x0.variable_like(2), Monomial::<3>::variable(2));
    }

    #[test]
    fn polynomial_helpers_match_receiver_ring() {
        let x0 = Poly::<3>::variable(0);

        assert_eq!(x0.zero_like(), Poly::<3>::zero());
        assert_eq!(x0.one_like(), Poly::<3>::one());
        assert_eq!(
            x0.constant_like(F::from(7)),
            Poly::<3>::constant(F::from(7))
        );
        assert_eq!(x0.variable_like(2), Poly::<3>::variable(2));
    }

    #[test]
    fn log2_ceil_is_the_smallest_exponent_that_fits_values() {
        assert_eq!(log2_ceil(0), 0);
        assert_eq!(log2_ceil(1), 0);
        assert_eq!(log2_ceil(2), 1);
        assert_eq!(log2_ceil(3), 2);
        assert_eq!(log2_ceil(4), 2);
        assert_eq!(log2_ceil(5), 3);
        assert_eq!(log2_ceil(8), 3);
        assert_eq!(log2_ceil(9), 4);
    }

    #[test]
    fn inner_product_multiplies_and_sums_pairs() {
        let lhs = [F::from(2), F::from(3), F::from(5)];
        let rhs = [F::from(7), F::from(11), F::from(13)];

        assert_eq!(inner_product(&lhs, &rhs), F::from(112));
    }
}
