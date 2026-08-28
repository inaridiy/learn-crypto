use super::helpers::lagrange_interpolation;
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;

pub trait UniPoly<F: Field> {
    fn degree(&self) -> usize;
    fn eval(&self, point: F) -> F;
}

/// `evals[i] = p(i)` という評価値で表した一変数多項式。
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct EvalsUniPoly<F: Field> {
    evals: Vec<F>,
}

impl<F: Field> EvalsUniPoly<F> {
    pub fn new(evals: Vec<F>) -> Self {
        assert!(
            !evals.is_empty(),
            "a polynomial needs at least one evaluation"
        );
        Self { evals }
    }

    pub fn evals(&self) -> &[F] {
        &self.evals
    }

    pub fn to_coeffs(&self) -> CoeffsUniPoly<F> {
        CoeffsUniPoly::from(self)
    }
}

impl<F: Field> UniPoly<F> for EvalsUniPoly<F> {
    fn degree(&self) -> usize {
        self.evals.len() - 1
    }

    fn eval(&self, point: F) -> F {
        self.to_coeffs().eval(point)
    }
}

/// `coeffs[i]` を `x^i` の係数として表した一変数多項式。
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize)]
pub struct CoeffsUniPoly<F: Field> {
    coeffs: Vec<F>,
}

impl<F: Field> CoeffsUniPoly<F> {
    pub fn new(coeffs: Vec<F>) -> Self {
        assert!(
            !coeffs.is_empty(),
            "a polynomial needs at least one coefficient"
        );
        Self { coeffs }
    }

    pub fn coeffs(&self) -> &[F] {
        &self.coeffs
    }
}

impl<F: Field> UniPoly<F> for CoeffsUniPoly<F> {
    fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    fn eval(&self, point: F) -> F {
        self.coeffs
            .iter()
            .rev()
            .fold(F::zero(), |value, &coeff| value * point + coeff)
    }
}

impl<F: Field> From<&EvalsUniPoly<F>> for CoeffsUniPoly<F> {
    fn from(poly: &EvalsUniPoly<F>) -> Self {
        Self::new(lagrange_interpolation(poly.evals()))
    }
}

impl<F: Field> From<EvalsUniPoly<F>> for CoeffsUniPoly<F> {
    fn from(poly: EvalsUniPoly<F>) -> Self {
        Self::from(&poly)
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;

    use super::{CoeffsUniPoly, EvalsUniPoly, UniPoly};

    fn f(value: u64) -> Fr {
        Fr::from(value)
    }

    #[test]
    fn evaluations_convert_to_coefficients() {
        // p(x) = 3x^2 + 2x + 5, evaluated at x = 0, 1, 2.
        let evals_poly = EvalsUniPoly::new(vec![f(5), f(10), f(21)]);
        let coeffs_poly = evals_poly.to_coeffs();

        assert_eq!(evals_poly.degree(), 2);
        assert_eq!(evals_poly.evals(), &[f(5), f(10), f(21)]);
        assert_eq!(coeffs_poly.degree(), 2);
        assert_eq!(coeffs_poly.coeffs(), &[f(5), f(2), f(3)]);
    }

    #[test]
    fn both_representations_evaluate_to_the_same_value() {
        let evals_poly = EvalsUniPoly::new(vec![f(5), f(10), f(21)]);
        let coeffs_poly = CoeffsUniPoly::new(vec![f(5), f(2), f(3)]);

        assert_eq!(evals_poly.eval(f(7)), f(166));
        assert_eq!(coeffs_poly.eval(f(7)), f(166));
    }

    #[test]
    fn constant_polynomials_work() {
        let evals_poly = EvalsUniPoly::new(vec![f(9)]);
        let coeffs_poly = CoeffsUniPoly::new(vec![f(9)]);

        assert_eq!(evals_poly.degree(), 0);
        assert_eq!(coeffs_poly.degree(), 0);
        assert_eq!(evals_poly.eval(f(123)), f(9));
        assert_eq!(coeffs_poly.eval(f(123)), f(9));
    }

    #[test]
    #[should_panic(expected = "a polynomial needs at least one evaluation")]
    fn empty_evaluations_are_rejected() {
        let _ = EvalsUniPoly::<Fr>::new(Vec::new());
    }

    #[test]
    #[should_panic(expected = "a polynomial needs at least one coefficient")]
    fn empty_coefficients_are_rejected() {
        let _ = CoeffsUniPoly::<Fr>::new(Vec::new());
    }
}
