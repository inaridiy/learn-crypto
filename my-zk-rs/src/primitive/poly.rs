use core::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_ff::Field;

/// `N` 変数多項式の単項式。
///
/// `exps[i]` が変数 `x_i` の指数を表す。たとえば `Monomial::<3>::new([2, 0, 1])`
/// は `x0^2 * x2` を表す。
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Monomial<const N: usize> {
    exps: [usize; N],
}

impl<const N: usize> Monomial<N> {
    /// 指数配列から単項式を作る。
    #[inline]
    pub fn new(exps: [usize; N]) -> Self {
        Self { exps }
    }

    /// 変数 `x_index` を表す単項式を作る。
    #[inline]
    pub fn variable(index: usize) -> Self {
        assert!(index < N, "variable index out of bounds");

        let mut exps = [0; N];
        exps[index] = 1;
        Self { exps }
    }

    /// 乗法単位元 `1` を表す単項式を作る。
    #[inline]
    pub fn one() -> Self {
        Self { exps: [0; N] }
    }

    /// 全変数の指数配列を返す。
    #[inline]
    pub fn exponents(&self) -> &[usize; N] {
        &self.exps
    }

    /// 変数 `x_index` の指数を返す。
    #[inline]
    pub fn exponent(&self, index: usize) -> usize {
        self.exps[index]
    }

    /// 単項式の全次数を返す。
    ///
    /// `x0^2 * x2` なら `3`。
    #[inline]
    pub fn degree(&self) -> usize {
        self.exps.iter().sum()
    }

    /// 単項式が `1` かどうかを返す。
    #[inline]
    pub fn is_one(&self) -> bool {
        self.exps.iter().all(|&x| x == 0)
    }

    /// 2つの単項式を掛ける。
    ///
    /// 同じ変数の指数を足し合わせる。
    #[inline]
    pub fn mul(&self, rhs: &Self) -> Self {
        let exps = std::array::from_fn(|i| self.exps[i] + rhs.exps[i]);
        Monomial { exps }
    }

    /// 点 `point = [x0, x1, ...]` で単項式を評価する。
    #[inline]
    pub fn eval<F: Field>(&self, point: &[F; N]) -> F {
        let mut value = F::one();
        for i in 0..N {
            match self.exps[i] {
                0 => {}
                1 => value *= point[i],
                exp => value *= point[i].pow([exp as u64]),
            }
        }

        value
    }
}

impl<const N: usize> fmt::Display for Monomial<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_one() {
            return write!(f, "1");
        }

        let mut first = true;
        for (i, exp) in self.exps.iter().enumerate() {
            if *exp == 0 {
                continue;
            }

            if !first {
                write!(f, "*")?;
            }

            if *exp == 1 {
                write!(f, "x{i}")?;
            } else {
                write!(f, "x{i}^{exp}")?;
            }
            first = false;
        }

        Ok(())
    }
}

/// `N` 変数多項式。
///
/// 内部表現は単項式順に正規化された `(Monomial, coeff)` の配列。係数が `0` に
/// なった項は保持しない。変数名は表示上 `x0, x1, ...` として扱う。
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct MvPolynomial<F: Field, const N: usize> {
    terms: Vec<(Monomial<N>, F)>,
}

impl<F: Field, const N: usize> MvPolynomial<F, N> {
    /// 零多項式を作る。
    #[inline]
    pub fn zero() -> Self {
        Self { terms: Vec::new() }
    }

    /// 定数多項式 `1` を作る。
    #[inline]
    pub fn one() -> Self {
        Self::constant(F::one())
    }

    /// 定数多項式を作る。
    ///
    /// `coeff == 0` の場合は零多項式になる。
    #[inline]
    pub fn constant(coeff: F) -> Self {
        Self::term(Monomial::one(), coeff)
    }

    /// 変数 `x_index` を表す多項式を作る。
    #[inline]
    pub fn variable(index: usize) -> Self {
        Self::term(Monomial::variable(index), F::one())
    }

    /// 1つの項だけを持つ多項式を作る。
    ///
    /// `coeff == 0` の場合は項を持たない零多項式になる。
    #[inline]
    pub fn term(monomial: Monomial<N>, coeff: F) -> Self {
        if coeff.is_zero() {
            Self::zero()
        } else {
            Self {
                terms: vec![(monomial, coeff)],
            }
        }
    }

    /// 項の列から多項式を作る。
    ///
    /// 同じ単項式の係数は足し合わされ、結果が `0` になった項は削除される。
    pub fn from_terms(terms: impl IntoIterator<Item = (Monomial<N>, F)>) -> Self {
        Self {
            terms: Self::normalize_terms(terms.into_iter().collect()),
        }
    }

    /// 保持している項を単項式順に走査する。
    #[inline]
    pub fn terms(&self) -> impl Iterator<Item = (&Monomial<N>, &F)> {
        self.terms.iter().map(|(monomial, coeff)| (monomial, coeff))
    }

    /// 係数が非零の項数を返す。
    #[inline]
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// 零多項式かどうかを返す。
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// 定数多項式 `1` かどうかを返す。
    #[inline]
    pub fn is_one_polynomial(&self) -> bool {
        self.terms.len() == 1 && self.terms[0].0.is_one() && self.terms[0].1 == F::one()
    }

    /// 多項式の全次数を返す。
    ///
    /// 各項の全次数 `sum_i exponent_i` の最大値。零多項式では `None`。
    /// SumCheck などで「各変数ごとの次数」が必要な場合は [`Self::degree_of`] または
    /// [`Self::degrees`] を使う。
    #[inline]
    pub fn degree(&self) -> Option<usize> {
        self.terms
            .iter()
            .map(|(monomial, _)| monomial.degree())
            .max()
    }

    /// 変数 `x_index` に関する次数を返す。
    ///
    /// 各項における `x_index` の指数の最大値。零多項式では `None`。
    /// たとえば `3*x0^2 + 2*x1 + 5` では `degree_of(0) == Some(2)`、
    /// `degree_of(1) == Some(1)`。
    #[inline]
    pub fn degree_of(&self, index: usize) -> Option<usize> {
        assert!(index < N, "variable index out of bounds");
        self.terms
            .iter()
            .map(|(monomial, _)| monomial.exponent(index))
            .max()
    }

    /// 各変数に関する次数をまとめて返す。
    ///
    /// 戻り値の `i` 番目が `x_i` に関する次数。零多項式では全要素が `0`。
    /// 零多項式と定数多項式を区別したい場合は [`Self::is_zero`] も確認する。
    #[inline]
    pub fn degrees(&self) -> [usize; N] {
        let mut degrees = [0; N];
        for (monomial, _) in &self.terms {
            for (i, exp) in monomial.exponents().iter().enumerate() {
                degrees[i] = degrees[i].max(*exp);
            }
        }
        degrees
    }

    /// multilinear polynomial かどうかを返す。
    ///
    /// すべての変数 `x_i` について `deg_i <= 1` のとき true。
    #[inline]
    pub fn is_multilinear(&self) -> bool {
        self.degrees().iter().all(|degree| *degree <= 1)
    }

    /// 項を追加する。
    ///
    /// 既に同じ単項式がある場合は係数を足す。合計が `0` になった項は削除する。
    #[inline]
    pub fn add_term(&mut self, monomial: Monomial<N>, coeff: F) {
        if coeff.is_zero() {
            return;
        }

        match self
            .terms
            .binary_search_by(|(current, _)| current.cmp(&monomial))
        {
            Ok(pos) => {
                self.terms[pos].1 += coeff;
                if self.terms[pos].1.is_zero() {
                    self.terms.remove(pos);
                }
            }
            Err(pos) => {
                self.terms.insert(pos, (monomial, coeff));
            }
        }
    }

    /// 点 `point = [x0, x1, ...]` で多項式を評価する。
    #[inline]
    pub fn eval(&self, point: &[F; N]) -> F {
        if self.is_zero() {
            return F::zero();
        }
        if self.terms.len() == 1 {
            let (monomial, coeff) = &self.terms[0];
            if monomial.is_one() {
                return *coeff;
            }
        }

        let degrees = self.degrees();
        let powers: Vec<Vec<F>> = (0..N)
            .map(|i| Self::scalar_powers(point[i], degrees[i]))
            .collect();

        let mut value = F::zero();
        for (monomial, coeff) in &self.terms {
            let mut term = *coeff;
            for (i, exp) in monomial.exponents().iter().enumerate() {
                if *exp != 0 {
                    term *= powers[i][*exp];
                }
            }
            value += term;
        }

        value
    }

    /// 多項式全体にスカラーを掛ける。
    ///
    /// `scalar == 0` の場合は零多項式を返す。
    #[inline]
    pub fn scale(mut self, scalar: F) -> Self {
        if scalar.is_zero() {
            return Self::zero();
        }
        if scalar == F::one() {
            return self;
        }

        for (_, coeff) in &mut self.terms {
            *coeff *= scalar;
        }
        self
    }

    /// 多項式の非負整数冪を計算する。
    #[inline]
    pub fn pow(&self, mut exp: usize) -> Self {
        if exp == 0 {
            return Self::one();
        }
        if exp == 1 {
            return self.clone();
        }
        if self.is_zero() {
            return Self::zero();
        }
        if self.is_one_polynomial() {
            return Self::one();
        }

        let mut base = self.clone();
        let mut result = Self::one();

        while exp > 0 {
            if exp & 1 == 1 {
                result *= &base;
            }

            exp >>= 1;
            if exp > 0 {
                base = &base * &base;
            }
        }

        result
    }

    /// 変数 `x_index` に値 `value` を代入し、その変数を消して `N - 1` 変数多項式にする。
    ///
    /// SumCheck の各ラウンドで verifier challenge を代入して次の多項式へ進む用途では、
    /// 通常このメソッドを使う。評価点の配列長も1つ短くなる。
    ///
    /// 例: `P(x0, x1, x2).curry_variable(1, r)` は
    /// `Q(x0, x1) = P(x0, r, x1)` を返す。
    #[inline]
    pub fn curry_variable(&self, index: usize, value: F) -> MvPolynomial<F, { N - 1 }>
    where
        [(); N - 1]:,
    {
        assert!(N > 0, "cannot drop a variable from a zero-var polynomial");
        assert!(index < N, "variable index out of bounds");

        let mut powers = vec![F::one()];
        let mut result = MvPolynomial::<F, { N - 1 }>::zero();
        for (monomial, coeff) in &self.terms {
            let mut exps = [0; N - 1];
            for old_i in 0..N {
                if old_i == index {
                    continue;
                }

                let new_i = if old_i < index { old_i } else { old_i - 1 };
                exps[new_i] = monomial.exponent(old_i);
            }

            let fixed_factor =
                Self::scalar_power_cached(&mut powers, value, monomial.exponent(index));
            result.add_term(Monomial::new(exps), *coeff * fixed_factor);
        }

        result
    }

    /// 先頭 `K` 個の変数を `values` で固定し、残りの `M` 変数多項式にする。
    ///
    /// 例: `P(x0, x1, x2, x3).curry_prefix(&[a, b])` は
    /// `Q(x0, x1) = P(a, b, x0, x1)` を返す。
    #[inline]
    pub fn curry_prefix<const K: usize, const M: usize>(
        &self,
        values: &[F; K],
    ) -> MvPolynomial<F, M> {
        assert_eq!(K + M, N, "invalid number of curried variables");

        let mut powers = vec![vec![F::one()]; K];
        let mut result = MvPolynomial::<F, M>::zero();
        for (monomial, coeff) in &self.terms {
            let mut fixed_factor = F::one();
            for i in 0..K {
                fixed_factor *=
                    Self::scalar_power_cached(&mut powers[i], values[i], monomial.exponent(i));
            }

            let mut exps = [0; M];
            for old_i in K..N {
                exps[old_i - K] = monomial.exponent(old_i);
            }

            result.add_term(Monomial::new(exps), *coeff * fixed_factor);
        }

        result
    }

    /// 末尾 `K` 個の変数を `values` で固定し、残りの `M` 変数多項式にする。
    ///
    /// 例: `P(x0, x1, x2, x3).curry_suffix(&[a, b])` は
    /// `Q(x0, x1) = P(x0, x1, a, b)` を返す。
    #[inline]
    pub fn curry_suffix<const K: usize, const M: usize>(
        &self,
        values: &[F; K],
    ) -> MvPolynomial<F, M> {
        assert_eq!(M + K, N, "invalid number of curried variables");

        let fixed_start = M;
        let mut powers = vec![vec![F::one()]; K];
        let mut result = MvPolynomial::<F, M>::zero();
        for (monomial, coeff) in &self.terms {
            let mut fixed_factor = F::one();
            for i in 0..K {
                let old_i = fixed_start + i;
                fixed_factor *=
                    Self::scalar_power_cached(&mut powers[i], values[i], monomial.exponent(old_i));
            }

            let mut exps = [0; M];
            for (old_i, exp) in exps.iter_mut().enumerate() {
                *exp = monomial.exponent(old_i);
            }

            result.add_term(Monomial::new(exps), *coeff * fixed_factor);
        }

        result
    }

    #[inline]
    fn scalar_powers(base: F, max_degree: usize) -> Vec<F> {
        let mut powers = Vec::with_capacity(max_degree + 1);
        powers.push(F::one());

        for exp in 1..=max_degree {
            powers.push(powers[exp - 1] * base);
        }

        powers
    }

    #[inline]
    fn scalar_power_cached(powers: &mut Vec<F>, base: F, exp: usize) -> F {
        while powers.len() <= exp {
            let next = *powers.last().expect("powers has one element") * base;
            powers.push(next);
        }

        powers[exp]
    }

    fn normalize_terms(mut terms: Vec<(Monomial<N>, F)>) -> Vec<(Monomial<N>, F)> {
        terms.retain(|(_, coeff)| !coeff.is_zero());
        terms.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0));

        let mut normalized: Vec<(Monomial<N>, F)> = Vec::with_capacity(terms.len());
        for (monomial, coeff) in terms {
            if let Some((last_monomial, last_coeff)) = normalized.last_mut() {
                if last_monomial == &monomial {
                    *last_coeff += coeff;
                    if last_coeff.is_zero() {
                        normalized.pop();
                    }
                    continue;
                }
            }

            normalized.push((monomial, coeff));
        }

        normalized
    }

    fn merge_terms(
        lhs: Vec<(Monomial<N>, F)>,
        rhs: Vec<(Monomial<N>, F)>,
        subtract_rhs: bool,
    ) -> Vec<(Monomial<N>, F)> {
        let mut merged = Vec::with_capacity(lhs.len() + rhs.len());
        let mut lhs = lhs.into_iter().peekable();
        let mut rhs = rhs.into_iter().peekable();

        loop {
            match (lhs.peek(), rhs.peek()) {
                (Some((lhs_monomial, _)), Some((rhs_monomial, _))) => {
                    match lhs_monomial.cmp(rhs_monomial) {
                        std::cmp::Ordering::Less => {
                            merged.push(lhs.next().expect("peeked"));
                        }
                        std::cmp::Ordering::Greater => {
                            let (monomial, coeff) = rhs.next().expect("peeked");
                            merged.push((monomial, if subtract_rhs { -coeff } else { coeff }));
                        }
                        std::cmp::Ordering::Equal => {
                            let (monomial, lhs_coeff) = lhs.next().expect("peeked");
                            let (_, rhs_coeff) = rhs.next().expect("peeked");
                            let coeff = if subtract_rhs {
                                lhs_coeff - rhs_coeff
                            } else {
                                lhs_coeff + rhs_coeff
                            };
                            if !coeff.is_zero() {
                                merged.push((monomial, coeff));
                            }
                        }
                    }
                }
                (Some(_), None) => {
                    merged.extend(lhs);
                    break;
                }
                (None, Some(_)) => {
                    if subtract_rhs {
                        merged.extend(rhs.map(|(monomial, coeff)| (monomial, -coeff)));
                    } else {
                        merged.extend(rhs);
                    }
                    break;
                }
                (None, None) => break,
            }
        }

        merged
    }
}

impl<F: Field + fmt::Display, const N: usize> fmt::Display for MvPolynomial<F, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let mut first = true;
        for (monomial, coeff) in &self.terms {
            if !first {
                write!(f, " + ")?;
            }

            if monomial.is_one() {
                write!(f, "{coeff}")?;
            } else if *coeff == F::one() {
                write!(f, "{monomial}")?;
            } else {
                write!(f, "{coeff}*{monomial}")?;
            }

            first = false;
        }

        Ok(())
    }
}

impl<F: Field, const N: usize> Default for MvPolynomial<F, N> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl<F: Field, const N: usize> AddAssign<&Self> for MvPolynomial<F, N> {
    #[inline]
    fn add_assign(&mut self, rhs: &Self) {
        self.terms = Self::merge_terms(std::mem::take(&mut self.terms), rhs.terms.clone(), false);
    }
}

impl<F: Field, const N: usize> AddAssign<Self> for MvPolynomial<F, N> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.terms = Self::merge_terms(std::mem::take(&mut self.terms), rhs.terms, false);
    }
}

impl<F: Field, const N: usize> Add<&Self> for MvPolynomial<F, N> {
    type Output = Self;

    #[inline]
    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<F: Field, const N: usize> Add<Self> for MvPolynomial<F, N> {
    type Output = Self;

    #[inline]
    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<F: Field, const N: usize> Add<&MvPolynomial<F, N>> for &MvPolynomial<F, N> {
    type Output = MvPolynomial<F, N>;

    #[inline]
    fn add(self, rhs: &MvPolynomial<F, N>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<F: Field, const N: usize> SubAssign<&Self> for MvPolynomial<F, N> {
    #[inline]
    fn sub_assign(&mut self, rhs: &Self) {
        self.terms = Self::merge_terms(std::mem::take(&mut self.terms), rhs.terms.clone(), true);
    }
}

impl<F: Field, const N: usize> SubAssign<Self> for MvPolynomial<F, N> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.terms = Self::merge_terms(std::mem::take(&mut self.terms), rhs.terms, true);
    }
}

impl<F: Field, const N: usize> Sub<&Self> for MvPolynomial<F, N> {
    type Output = Self;

    #[inline]
    fn sub(mut self, rhs: &Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<F: Field, const N: usize> Sub<Self> for MvPolynomial<F, N> {
    type Output = Self;

    #[inline]
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<F: Field, const N: usize> Sub<&MvPolynomial<F, N>> for &MvPolynomial<F, N> {
    type Output = MvPolynomial<F, N>;

    #[inline]
    fn sub(self, rhs: &MvPolynomial<F, N>) -> Self::Output {
        self.clone() - rhs
    }
}

impl<F: Field, const N: usize> Neg for MvPolynomial<F, N> {
    type Output = Self;

    #[inline]
    fn neg(mut self) -> Self::Output {
        for (_, coeff) in &mut self.terms {
            *coeff = -*coeff;
        }
        self
    }
}

impl<F: Field, const N: usize> Neg for &MvPolynomial<F, N> {
    type Output = MvPolynomial<F, N>;

    #[inline]
    fn neg(self) -> Self::Output {
        self.clone().neg()
    }
}

impl<F: Field, const N: usize> MulAssign<&Self> for MvPolynomial<F, N> {
    #[inline]
    fn mul_assign(&mut self, rhs: &Self) {
        *self = &*self * rhs;
    }
}

impl<F: Field, const N: usize> MulAssign<Self> for MvPolynomial<F, N> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl<F: Field, const N: usize> Mul<&Self> for MvPolynomial<F, N> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: &Self) -> Self::Output {
        &self * rhs
    }
}

impl<F: Field, const N: usize> Mul<Self> for MvPolynomial<F, N> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        if self.is_one_polynomial() {
            return rhs;
        }
        if rhs.is_one_polynomial() {
            return self;
        }

        &self * &rhs
    }
}

impl<F: Field, const N: usize> Mul<&MvPolynomial<F, N>> for &MvPolynomial<F, N> {
    type Output = MvPolynomial<F, N>;

    #[inline]
    fn mul(self, rhs: &MvPolynomial<F, N>) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return MvPolynomial::zero();
        }
        if self.is_one_polynomial() {
            return rhs.clone();
        }
        if rhs.is_one_polynomial() {
            return self.clone();
        }

        let mut terms = Vec::with_capacity(self.num_terms().saturating_mul(rhs.num_terms()));
        for (lhs_monomial, lhs_coeff) in &self.terms {
            for (rhs_monomial, rhs_coeff) in &rhs.terms {
                let coeff = *lhs_coeff * *rhs_coeff;
                if !coeff.is_zero() {
                    terms.push((lhs_monomial.mul(rhs_monomial), coeff));
                }
            }
        }

        terms.sort_unstable_by(|lhs, rhs| lhs.0.cmp(&rhs.0));

        let mut product_terms = Vec::with_capacity(terms.len());
        let mut terms = terms.into_iter();
        if let Some((mut current_monomial, mut current_coeff)) = terms.next() {
            for (monomial, coeff) in terms {
                if monomial == current_monomial {
                    current_coeff += coeff;
                } else {
                    if !current_coeff.is_zero() {
                        product_terms.push((current_monomial, current_coeff));
                    }
                    current_monomial = monomial;
                    current_coeff = coeff;
                }
            }

            if !current_coeff.is_zero() {
                product_terms.push((current_monomial, current_coeff));
            }
        }

        MvPolynomial {
            terms: product_terms,
        }
    }
}

impl<F: Field, const N: usize> Mul<F> for MvPolynomial<F, N> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: F) -> Self::Output {
        self.scale(rhs)
    }
}

impl<F: Field, const N: usize> Mul<F> for &MvPolynomial<F, N> {
    type Output = MvPolynomial<F, N>;

    #[inline]
    fn mul(self, rhs: F) -> Self::Output {
        self.clone().scale(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::{Monomial, MvPolynomial};
    use ark_bls12_381::Fr as F;

    type Poly<const N: usize> = MvPolynomial<F, N>;

    fn f(x: u64) -> F {
        F::from(x)
    }

    #[test]
    fn monomial_evaluates_at_point() {
        let monomial = Monomial::<2>::new([2, 1]);

        assert_eq!(monomial.eval(&[f(3), f(4)]), f(36));
    }

    #[test]
    fn construction_combines_like_terms_and_removes_zero_coefficients() {
        let x = Monomial::<2>::variable(0);
        let poly = Poly::from_terms([(x.clone(), f(3)), (x, -f(3)), (Monomial::one(), f(5))]);

        assert_eq!(poly.num_terms(), 1);
        assert_eq!(poly.eval(&[f(10), f(20)]), f(5));
    }

    #[test]
    fn polynomial_evaluates_at_point() {
        let poly = Poly::from_terms([
            (Monomial::new([2, 0]), f(3)),
            (Monomial::new([0, 1]), f(2)),
            (Monomial::one(), f(5)),
        ]);

        assert_eq!(poly.degree(), Some(2));
        assert_eq!(poly.degree_of(0), Some(2));
        assert_eq!(poly.degree_of(1), Some(1));
        assert_eq!(poly.degrees(), [2, 1]);
        assert_eq!(poly.eval(&[f(2), f(7)]), f(31));
    }

    #[test]
    fn addition_subtraction_and_multiplication_work() {
        let x = Poly::<2>::variable(0);
        let y = Poly::<2>::variable(1);
        let one = Poly::<2>::one();

        let product = (&x + &one) * (&y + &one);

        assert_eq!(product.num_terms(), 4);
        assert_eq!(product.degree(), Some(2));
        assert_eq!(product.eval(&[f(2), f(3)]), f(12));
        assert!((product.clone() - product).is_zero());
    }

    #[test]
    fn polynomial_power_uses_repeated_multiplication() {
        let x = Poly::<2>::variable(0);
        let one = Poly::<2>::one();

        let poly = (&x + &one).pow(3);

        assert_eq!(poly.eval(&[f(2), f(99)]), f(27));
    }

    #[test]
    fn curry_variable_evaluates_one_variable_and_removes_it() {
        let poly = Poly::from_terms([
            (Monomial::new([2, 0]), f(3)),
            (Monomial::new([0, 1]), f(2)),
            (Monomial::one(), f(5)),
        ]);

        let pinned_x = poly.curry_variable(0, f(2));
        let pinned_y = poly.curry_variable(1, f(7));

        assert_eq!(pinned_x.degrees(), [1]);
        assert_eq!(pinned_y.degrees(), [2]);
        assert_eq!(pinned_x.eval(&[f(7)]), f(31));
        assert_eq!(pinned_y.eval(&[f(2)]), f(31));
    }

    #[test]
    fn curry_prefix_evaluates_leading_variables_and_removes_them() {
        let poly = Poly::from_terms([
            (Monomial::new([1, 0, 1, 0]), f(2)),
            (Monomial::new([0, 1, 0, 1]), f(3)),
            (Monomial::one(), f(5)),
        ]);

        let curried = poly.curry_prefix(&[f(7), f(11)]);

        assert_eq!(curried.degrees(), [1, 1]);
        assert_eq!(
            curried.eval(&[f(13), f(17)]),
            poly.eval(&[f(7), f(11), f(13), f(17)])
        );
    }

    #[test]
    fn curry_suffix_evaluates_trailing_variables_and_removes_them() {
        let poly = Poly::from_terms([
            (Monomial::new([1, 0, 1, 0]), f(2)),
            (Monomial::new([0, 1, 0, 1]), f(3)),
            (Monomial::one(), f(5)),
        ]);

        let curried = poly.curry_suffix(&[f(13), f(17)]);

        assert_eq!(curried.degrees(), [1, 1]);
        assert_eq!(
            curried.eval(&[f(7), f(11)]),
            poly.eval(&[f(7), f(11), f(13), f(17)])
        );
    }
}
