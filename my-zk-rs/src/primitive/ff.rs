use std::{
    cmp::Ordering,
    fmt::{Debug, Display, Formatter, Result as FmtResult},
    hash::{Hash, Hasher},
    iter::{Product, Sum},
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
};

use ark_ff::{
    AdditiveGroup, BigInt as ArkBigInt, FftField, Field, LegendreSymbol, One, PrimeField,
    SqrtPrecomputation, Zero,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, Compress, EmptyFlags, Flags, SerializationError, Valid, Validate,
    buffer_byte_size,
};
use ark_std::rand::{
    Rng,
    distributions::{Distribution, Standard},
};
use num_bigint::BigUint;
use num_traits::Signed;
use zeroize::Zeroize;

/// Configuration for a Montgomery prime field over `N` little-endian `u64` limbs.
pub trait FpConfig<const N: usize>: Send + Sync + 'static + Sized {
    /// Little-endian limbs of the prime modulus.
    const MODULUS: [u64; N];

    /// A multiplicative generator, represented in canonical little-endian limbs.
    const GENERATOR: [u64; N];
}

#[derive(Debug)]
pub struct F25519Config;

impl FpConfig<4> for F25519Config {
    // P = 2^255 - 19
    const MODULUS: [u64; 4] = [
        0xffffffffffffffed,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x7fffffffffffffff,
    ];

    const GENERATOR: [u64; 4] = [2, 0, 0, 0];
}

pub type Fp25519 = Fp<F25519Config, 4>;

/// Number of 64-bit limbs used by [`Fp25519`].
pub const NUM_LIMBS: usize = 4;

/// Little-endian limbs of the modulus `2^255 - 19`.
pub const MODULUS_LIMBS: [u64; NUM_LIMBS] = F25519Config::MODULUS;

pub struct Fp<C: FpConfig<N>, const N: usize> {
    limbs: [u64; N],
    _config: PhantomData<C>,
}

impl<C: FpConfig<N>, const N: usize> Copy for Fp<C, N> {}

impl<C: FpConfig<N>, const N: usize> Clone for Fp<C, N> {
    #[inline(always)]
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: FpConfig<N>, const N: usize> Fp<C, N> {
    pub const MODULUS: [u64; N] = C::MODULUS;
    pub const GENERATOR_LIMBS: [u64; N] = C::GENERATOR;
    pub const MONT_INV: u64 = mont_inv(C::MODULUS);
    pub const R: [u64; N] = mont_r(C::MODULUS);
    pub const R2: [u64; N] = mont_r2(C::MODULUS);
    pub const TWO_ADICITY: u32 = ArkBigInt(C::MODULUS).two_adic_valuation();
    pub const TRACE: ArkBigInt<N> = ArkBigInt(C::MODULUS).two_adic_coefficient();
    pub const TWO_ADIC_ROOT_OF_UNITY_LIMBS: [u64; N] = to_montgomery(
        const_pow_mod(C::GENERATOR, Self::TRACE.0, C::MODULUS),
        C::MODULUS,
    );
    pub const ZERO: Self = Self::from_raw_montgomery([0u64; N]);
    pub const ONE: Self = Self::from_raw_montgomery(Self::R);
    pub const NEG_ONE: Self = Self::from_raw_montgomery(mont_neg_one(C::MODULUS));
    pub const GENERATOR: Self = Self::from_raw_montgomery(to_montgomery(C::GENERATOR, C::MODULUS));
    pub const TWO_ADIC_ROOT_OF_UNITY: Self =
        Self::from_raw_montgomery(Self::TWO_ADIC_ROOT_OF_UNITY_LIMBS);

    #[inline(always)]
    pub const fn from_raw_montgomery(limbs: [u64; N]) -> Self {
        Self {
            limbs,
            _config: PhantomData,
        }
    }

    #[inline(always)]
    pub fn new(x: [u64; N]) -> Self {
        let reduced = reduce_limbs(x, C::MODULUS);
        Self::from_raw_montgomery(mont_mul(reduced, Self::R2, C::MODULUS, Self::MONT_INV))
    }

    #[inline(always)]
    pub fn zero() -> Self {
        Self::ZERO
    }

    #[inline(always)]
    pub fn one() -> Self {
        Self::ONE
    }

    #[inline(always)]
    pub fn from_u64(x: u64) -> Self {
        if x == 0 {
            return Self::ZERO;
        }
        if x == 1 {
            return Self::ONE;
        }

        let mut limbs = [0u64; N];
        limbs[0] = x;
        let reduced = reduce_limbs(limbs, C::MODULUS);
        if high_limbs_are_zero(&reduced) {
            Self::from_raw_montgomery(mont_mul_limb(
                reduced[0],
                Self::R2,
                C::MODULUS,
                Self::MONT_INV,
            ))
        } else {
            Self::from_raw_montgomery(mont_mul(reduced, Self::R2, C::MODULUS, Self::MONT_INV))
        }
    }

    #[inline(always)]
    pub fn to_limbs(self) -> [u64; N] {
        mont_reduce(self.limbs, C::MODULUS, Self::MONT_INV)
    }

    #[inline(always)]
    pub fn add(self, rhs: Self) -> Self {
        let (mut sum, carry) = add_raw(self.limbs, rhs.limbs);
        if carry || cmp_limbs(&sum, &C::MODULUS) != Ordering::Less {
            sum = sub_raw(sum, C::MODULUS).0;
        }

        Self::from_raw_montgomery(sum)
    }

    #[inline(always)]
    pub fn sub(self, rhs: Self) -> Self {
        let (mut diff, borrow) = sub_raw(self.limbs, rhs.limbs);
        if borrow {
            diff = add_raw(diff, C::MODULUS).0;
        }

        Self::from_raw_montgomery(diff)
    }

    #[inline(always)]
    pub fn mul(self, rhs: Self) -> Self {
        Self::from_raw_montgomery(mont_mul(self.limbs, rhs.limbs, C::MODULUS, Self::MONT_INV))
    }

    #[inline(always)]
    pub fn square(self) -> Self {
        self.mul(self)
    }

    #[inline(always)]
    pub fn double(self) -> Self {
        self.add(self)
    }

    #[inline(always)]
    pub fn inverse(self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(self.pow(modulus_minus_two(C::MODULUS)))
        }
    }

    #[inline]
    pub fn pow<S: AsRef<[u64]>>(self, exp: S) -> Self {
        let mut result = Self::one();
        let mut base = self;

        for &limb in exp.as_ref() {
            let mut bits = limb;
            for _ in 0..64 {
                if bits & 1 == 1 {
                    result *= base;
                }
                base.square_in_place();
                bits >>= 1;
            }
        }

        result
    }

    #[inline(always)]
    pub fn mul_batch(lhs: &[Self], rhs: &[Self], out: &mut [Self]) {
        assert_eq!(lhs.len(), rhs.len());
        assert_eq!(lhs.len(), out.len());

        for ((out, lhs), rhs) in out.iter_mut().zip(lhs).zip(rhs) {
            *out = (*lhs).mul(*rhs);
        }
    }

    #[inline(always)]
    fn is_zero_montgomery(&self) -> bool {
        self.limbs.iter().all(|&limb| limb == 0)
    }

    #[inline]
    fn sqrt_tonelli_shanks(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero());
        }
        if self.legendre() != LegendreSymbol::QuadraticResidue {
            return None;
        }

        let q = Self::TRACE.0;
        let mut q_plus_one_div_two = q;
        let (sum, carry) = add_raw(q_plus_one_div_two, one_limbs());
        debug_assert!(!carry);
        q_plus_one_div_two = div2(sum);

        let mut c = Self::GENERATOR.pow(q);
        let mut x = (*self).pow(q_plus_one_div_two);
        let mut t = (*self).pow(q);
        let mut m = Self::TWO_ADICITY;

        while !t.is_one() {
            let mut i = 1u32;
            let mut t2i = t.square();
            while !t2i.is_one() {
                t2i.square_in_place();
                i += 1;
                if i >= m {
                    return None;
                }
            }

            let mut b = c;
            for _ in 0..(m - i - 1) {
                b.square_in_place();
            }

            x *= b;
            c = b.square();
            t *= c;
            m = i;
        }

        (x.square() == *self).then_some(x)
    }
}

impl<C: FpConfig<N>, const N: usize> Default for Fp<C, N> {
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<C: FpConfig<N>, const N: usize> PartialEq for Fp<C, N> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.limbs == other.limbs
    }
}

impl<C: FpConfig<N>, const N: usize> Eq for Fp<C, N> {}

impl<C: FpConfig<N>, const N: usize> Hash for Fp<C, N> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.to_limbs().hash(state);
    }
}

impl<C: FpConfig<N>, const N: usize> Debug for Fp<C, N> {
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        Debug::fmt(&self.into_bigint(), f)
    }
}

impl<C: FpConfig<N>, const N: usize> Display for Fp<C, N> {
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        Display::fmt(&self.into_bigint(), f)
    }
}

impl<C: FpConfig<N>, const N: usize> Ord for Fp<C, N> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        cmp_limbs(&self.to_limbs(), &other.to_limbs())
    }
}

impl<C: FpConfig<N>, const N: usize> PartialOrd for Fp<C, N> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<C: FpConfig<N>, const N: usize> Zero for Fp<C, N> {
    #[inline(always)]
    fn zero() -> Self {
        Self::zero()
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.is_zero_montgomery()
    }
}

impl<C: FpConfig<N>, const N: usize> One for Fp<C, N> {
    #[inline(always)]
    fn one() -> Self {
        Self::one()
    }

    #[inline(always)]
    fn is_one(&self) -> bool {
        self.limbs == Self::R
    }
}

impl<C: FpConfig<N>, const N: usize> AdditiveGroup for Fp<C, N> {
    type Scalar = Self;

    const ZERO: Self = Self::ZERO;

    #[inline(always)]
    fn double(&self) -> Self {
        Fp::double(*self)
    }

    #[inline(always)]
    fn double_in_place(&mut self) -> &mut Self {
        *self = Fp::double(*self);
        self
    }

    #[inline(always)]
    fn neg_in_place(&mut self) -> &mut Self {
        *self = -*self;
        self
    }
}

impl<C: FpConfig<N>, const N: usize> Field for Fp<C, N> {
    type BasePrimeField = Self;

    const SQRT_PRECOMP: Option<SqrtPrecomputation<Self>> = None;
    const ONE: Self = Self::ONE;
    const NEG_ONE: Self = Self::NEG_ONE;

    #[inline(always)]
    fn characteristic() -> &'static [u64] {
        &C::MODULUS
    }

    #[inline(always)]
    fn extension_degree() -> u64 {
        1
    }

    #[inline(always)]
    fn to_base_prime_field_elements(&self) -> impl Iterator<Item = Self::BasePrimeField> {
        std::iter::once(*self)
    }

    #[inline]
    fn from_base_prime_field_elems(
        elems: impl IntoIterator<Item = Self::BasePrimeField>,
    ) -> Option<Self> {
        let mut elems = elems.into_iter();
        let first = elems.next()?;
        elems.next().is_none().then_some(first)
    }

    #[inline(always)]
    fn from_base_prime_field(elem: Self::BasePrimeField) -> Self {
        elem
    }

    #[inline]
    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        if F::BIT_SIZE > 8 {
            return None;
        }

        let output_byte_size = buffer_byte_size(Self::MODULUS_BIT_SIZE as usize + F::BIT_SIZE);
        let modulus_byte_size = buffer_byte_size(Self::MODULUS_BIT_SIZE as usize);
        let buffer_len = output_byte_size.max(modulus_byte_size);
        debug_assert!(buffer_len <= N * 8 + 1);

        let mut flags = F::default();
        let flags_pos = output_byte_size.saturating_sub(1);
        let mut limbs = [0u64; N];

        for pos in 0..buffer_len {
            let mut byte = bytes.get(pos).copied().unwrap_or(0);
            if output_byte_size != 0 && pos == flags_pos {
                flags = F::from_u8_remove_flags(&mut byte)?;
            }
            if pos < modulus_byte_size {
                limbs[pos / 8] |= (byte as u64) << ((pos % 8) * 8);
            }
        }

        mask_unused_limb_bits(&mut limbs, Self::MODULUS_BIT_SIZE as usize);
        Self::from_bigint(ArkBigInt(limbs)).map(|value| (value, flags))
    }

    #[inline]
    fn legendre(&self) -> LegendreSymbol {
        let s = (*self).pow(Self::MODULUS_MINUS_ONE_DIV_TWO);
        if s.is_zero() {
            LegendreSymbol::Zero
        } else if s.is_one() {
            LegendreSymbol::QuadraticResidue
        } else {
            LegendreSymbol::QuadraticNonResidue
        }
    }

    #[inline]
    fn sqrt(&self) -> Option<Self> {
        self.sqrt_tonelli_shanks()
    }

    #[inline(always)]
    fn square(&self) -> Self {
        Fp::square(*self)
    }

    #[inline(always)]
    fn square_in_place(&mut self) -> &mut Self {
        *self = Fp::square(*self);
        self
    }

    #[inline(always)]
    fn inverse(&self) -> Option<Self> {
        Fp::inverse(*self)
    }

    #[inline(always)]
    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        self.inverse().map(|inverse| {
            *self = inverse;
            self
        })
    }

    #[inline(always)]
    fn frobenius_map_in_place(&mut self, _: usize) {}

    #[inline(always)]
    fn mul_by_base_prime_field(&self, elem: &Self::BasePrimeField) -> Self {
        *self * elem
    }
}

impl<C: FpConfig<N>, const N: usize> FftField for Fp<C, N> {
    const GENERATOR: Self = Self::GENERATOR;
    const TWO_ADICITY: u32 = Self::TWO_ADICITY;
    const TWO_ADIC_ROOT_OF_UNITY: Self = Self::TWO_ADIC_ROOT_OF_UNITY;
}

impl<C: FpConfig<N>, const N: usize> PrimeField for Fp<C, N> {
    type BigInt = ArkBigInt<N>;

    const MODULUS: Self::BigInt = ArkBigInt(C::MODULUS);
    const MODULUS_MINUS_ONE_DIV_TWO: Self::BigInt = ArkBigInt(C::MODULUS).divide_by_2_round_down();
    const MODULUS_BIT_SIZE: u32 = ArkBigInt(C::MODULUS).const_num_bits();
    const TRACE: Self::BigInt = ArkBigInt(C::MODULUS).two_adic_coefficient();
    const TRACE_MINUS_ONE_DIV_TWO: Self::BigInt = ArkBigInt(C::MODULUS)
        .two_adic_coefficient()
        .divide_by_2_round_down();

    #[inline(always)]
    fn from_bigint(repr: Self::BigInt) -> Option<Self> {
        (repr < ArkBigInt(C::MODULUS)).then(|| {
            Self::from_raw_montgomery(mont_mul(repr.0, Self::R2, C::MODULUS, Self::MONT_INV))
        })
    }

    #[inline(always)]
    fn into_bigint(self) -> Self::BigInt {
        ArkBigInt(self.to_limbs())
    }
}

impl<C: FpConfig<N>, const N: usize> CanonicalSerializeWithFlags for Fp<C, N> {
    fn serialize_with_flags<W: ark_std::io::Write, F: Flags>(
        &self,
        mut writer: W,
        flags: F,
    ) -> Result<(), SerializationError> {
        if F::BIT_SIZE > 8 {
            return Err(SerializationError::NotEnoughSpace);
        }

        let output_byte_size = self.serialized_size_with_flags::<F>();
        debug_assert!(output_byte_size <= N * 8 + 1);
        let limbs = self.to_limbs();

        for pos in 0..output_byte_size {
            let mut byte = if pos < N * 8 {
                (limbs[pos / 8] >> ((pos % 8) * 8)) as u8
            } else {
                0
            };
            if pos + 1 == output_byte_size {
                byte |= flags.u8_bitmask();
            }
            writer.write_all(&[byte])?;
        }
        Ok(())
    }

    #[inline]
    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        buffer_byte_size(Self::MODULUS_BIT_SIZE as usize + F::BIT_SIZE)
    }
}

impl<C: FpConfig<N>, const N: usize> CanonicalSerialize for Fp<C, N> {
    #[inline]
    fn serialize_with_mode<W: ark_std::io::Write>(
        &self,
        writer: W,
        _compress: Compress,
    ) -> Result<(), SerializationError> {
        self.serialize_with_flags(writer, EmptyFlags)
    }

    #[inline]
    fn serialized_size(&self, _compress: Compress) -> usize {
        self.serialized_size_with_flags::<EmptyFlags>()
    }
}

impl<C: FpConfig<N>, const N: usize> CanonicalDeserializeWithFlags for Fp<C, N> {
    fn deserialize_with_flags<R: ark_std::io::Read, F: Flags>(
        mut reader: R,
    ) -> Result<(Self, F), SerializationError> {
        if F::BIT_SIZE > 8 {
            return Err(SerializationError::NotEnoughSpace);
        }

        let output_byte_size = buffer_byte_size(Self::MODULUS_BIT_SIZE as usize + F::BIT_SIZE);
        let modulus_byte_size = buffer_byte_size(Self::MODULUS_BIT_SIZE as usize);
        debug_assert!(output_byte_size <= N * 8 + 1);

        let mut flags = F::default();
        let mut limbs = [0u64; N];

        for pos in 0..output_byte_size {
            let mut byte = [0u8; 1];
            reader.read_exact(&mut byte)?;
            let mut byte = byte[0];
            if pos + 1 == output_byte_size {
                flags = F::from_u8_remove_flags(&mut byte)
                    .ok_or(SerializationError::UnexpectedFlags)?;
            }
            if pos < modulus_byte_size {
                limbs[pos / 8] |= (byte as u64) << ((pos % 8) * 8);
            }
        }

        mask_unused_limb_bits(&mut limbs, Self::MODULUS_BIT_SIZE as usize);
        Self::from_bigint(ArkBigInt(limbs))
            .map(|value| (value, flags))
            .ok_or(SerializationError::InvalidData)
    }
}

impl<C: FpConfig<N>, const N: usize> Valid for Fp<C, N> {
    const TRIVIAL_CHECK: bool = true;

    #[inline(always)]
    fn check(&self) -> Result<(), SerializationError> {
        Ok(())
    }
}

impl<C: FpConfig<N>, const N: usize> CanonicalDeserialize for Fp<C, N> {
    #[inline]
    fn deserialize_with_mode<R: ark_std::io::Read>(
        reader: R,
        _compress: Compress,
        _validate: Validate,
    ) -> Result<Self, SerializationError> {
        Self::deserialize_with_flags::<R, EmptyFlags>(reader).map(|(value, _)| value)
    }
}

impl<C: FpConfig<N>, const N: usize> Distribution<Fp<C, N>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Fp<C, N> {
        loop {
            let mut limbs = [0u64; N];
            for limb in &mut limbs {
                *limb = rng.next_u64();
            }

            let shave_bits = 64 * N - Fp::<C, N>::MODULUS_BIT_SIZE as usize;
            if shave_bits > 0 {
                limbs[N - 1] &= u64::MAX >> shave_bits;
            }

            if let Some(value) = Fp::<C, N>::from_bigint(ArkBigInt(limbs)) {
                return value;
            }
        }
    }
}

impl<C: FpConfig<N>, const N: usize> Zeroize for Fp<C, N> {
    #[inline]
    fn zeroize(&mut self) {
        self.limbs.zeroize();
    }
}

impl<C: FpConfig<N>, const N: usize> FromStr for Fp<C, N> {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let modulus = num_bigint::BigInt::from(ArkBigInt(C::MODULUS));
        let mut value = num_bigint::BigInt::from_str(s).map_err(|_| ())? % &modulus;
        if value.is_negative() {
            value += modulus;
        }

        let value = BigUint::try_from(value).map_err(|_| ())?;
        let value = ArkBigInt::<N>::try_from(value).map_err(|_| ())?;
        Self::from_bigint(value).ok_or(())
    }
}

impl<C: FpConfig<N>, const N: usize> From<ArkBigInt<N>> for Fp<C, N> {
    #[inline(always)]
    fn from(value: ArkBigInt<N>) -> Self {
        Self::from_bigint(value).expect("integer must be smaller than the field modulus")
    }
}

impl<C: FpConfig<N>, const N: usize> From<Fp<C, N>> for ArkBigInt<N> {
    #[inline(always)]
    fn from(value: Fp<C, N>) -> Self {
        value.into_bigint()
    }
}

impl<C: FpConfig<N>, const N: usize> From<BigUint> for Fp<C, N> {
    #[inline]
    fn from(value: BigUint) -> Self {
        Self::from_le_bytes_mod_order(&value.to_bytes_le())
    }
}

impl<C: FpConfig<N>, const N: usize> From<Fp<C, N>> for BigUint {
    #[inline]
    fn from(value: Fp<C, N>) -> Self {
        value.into_bigint().into()
    }
}

impl<C: FpConfig<N>, const N: usize> From<u128> for Fp<C, N> {
    #[inline]
    fn from(value: u128) -> Self {
        let mut limbs = [0u64; N];
        limbs[0] = value as u64;
        if N > 1 {
            limbs[1] = (value >> 64) as u64;
        }
        Self::new(limbs)
    }
}

impl<C: FpConfig<N>, const N: usize> From<i128> for Fp<C, N> {
    #[inline]
    fn from(value: i128) -> Self {
        let abs = Self::from(value.unsigned_abs());
        if value.is_negative() { -abs } else { abs }
    }
}

impl<C: FpConfig<N>, const N: usize> From<u64> for Fp<C, N> {
    #[inline(always)]
    fn from(value: u64) -> Self {
        Self::from_u64(value)
    }
}

impl<C: FpConfig<N>, const N: usize> From<i64> for Fp<C, N> {
    #[inline]
    fn from(value: i64) -> Self {
        let abs = Self::from(value.unsigned_abs());
        if value.is_negative() { -abs } else { abs }
    }
}

impl<C: FpConfig<N>, const N: usize> From<u32> for Fp<C, N> {
    #[inline(always)]
    fn from(value: u32) -> Self {
        Self::from(value as u64)
    }
}

impl<C: FpConfig<N>, const N: usize> From<i32> for Fp<C, N> {
    #[inline]
    fn from(value: i32) -> Self {
        let abs = Self::from(value.unsigned_abs());
        if value.is_negative() { -abs } else { abs }
    }
}

impl<C: FpConfig<N>, const N: usize> From<u16> for Fp<C, N> {
    #[inline(always)]
    fn from(value: u16) -> Self {
        Self::from(value as u64)
    }
}

impl<C: FpConfig<N>, const N: usize> From<i16> for Fp<C, N> {
    #[inline]
    fn from(value: i16) -> Self {
        let abs = Self::from(value.unsigned_abs());
        if value.is_negative() { -abs } else { abs }
    }
}

impl<C: FpConfig<N>, const N: usize> From<u8> for Fp<C, N> {
    #[inline(always)]
    fn from(value: u8) -> Self {
        Self::from(value as u64)
    }
}

impl<C: FpConfig<N>, const N: usize> From<i8> for Fp<C, N> {
    #[inline]
    fn from(value: i8) -> Self {
        let abs = Self::from(value.unsigned_abs());
        if value.is_negative() { -abs } else { abs }
    }
}

impl<C: FpConfig<N>, const N: usize> From<bool> for Fp<C, N> {
    #[inline(always)]
    fn from(value: bool) -> Self {
        Self::from(value as u64)
    }
}

impl<C: FpConfig<N>, const N: usize> Neg for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn neg(self) -> Self::Output {
        if self.is_zero() {
            self
        } else {
            Self::from_raw_montgomery(sub_raw(C::MODULUS, self.limbs).0)
        }
    }
}

impl<C: FpConfig<N>, const N: usize> Add<&Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn add(self, rhs: &Self) -> Self::Output {
        Fp::add(self, *rhs)
    }
}

impl<C: FpConfig<N>, const N: usize> Add<Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn add(self, rhs: Self) -> Self::Output {
        Fp::add(self, rhs)
    }
}

impl<C: FpConfig<N>, const N: usize> Add<&mut Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn add(self, rhs: &mut Self) -> Self::Output {
        self + &*rhs
    }
}

impl<C: FpConfig<N>, const N: usize> AddAssign<&Self> for Fp<C, N> {
    #[inline(always)]
    fn add_assign(&mut self, rhs: &Self) {
        *self = Fp::add(*self, *rhs);
    }
}

impl<C: FpConfig<N>, const N: usize> AddAssign<Self> for Fp<C, N> {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl<C: FpConfig<N>, const N: usize> AddAssign<&mut Self> for Fp<C, N> {
    #[inline(always)]
    fn add_assign(&mut self, rhs: &mut Self) {
        *self += &*rhs;
    }
}

impl<C: FpConfig<N>, const N: usize> Sub<&Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn sub(self, rhs: &Self) -> Self::Output {
        Fp::sub(self, *rhs)
    }
}

impl<C: FpConfig<N>, const N: usize> Sub<Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn sub(self, rhs: Self) -> Self::Output {
        Fp::sub(self, rhs)
    }
}

impl<C: FpConfig<N>, const N: usize> Sub<&mut Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn sub(self, rhs: &mut Self) -> Self::Output {
        self - &*rhs
    }
}

impl<C: FpConfig<N>, const N: usize> SubAssign<&Self> for Fp<C, N> {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: &Self) {
        *self = Fp::sub(*self, *rhs);
    }
}

impl<C: FpConfig<N>, const N: usize> SubAssign<Self> for Fp<C, N> {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs;
    }
}

impl<C: FpConfig<N>, const N: usize> SubAssign<&mut Self> for Fp<C, N> {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: &mut Self) {
        *self -= &*rhs;
    }
}

impl<C: FpConfig<N>, const N: usize> Mul<&Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn mul(self, rhs: &Self) -> Self::Output {
        Fp::mul(self, *rhs)
    }
}

impl<C: FpConfig<N>, const N: usize> Mul<Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn mul(self, rhs: Self) -> Self::Output {
        Fp::mul(self, rhs)
    }
}

impl<C: FpConfig<N>, const N: usize> Mul<&mut Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn mul(self, rhs: &mut Self) -> Self::Output {
        self * &*rhs
    }
}

impl<C: FpConfig<N>, const N: usize> MulAssign<&Self> for Fp<C, N> {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: &Self) {
        *self = Fp::mul(*self, *rhs);
    }
}

impl<C: FpConfig<N>, const N: usize> MulAssign<Self> for Fp<C, N> {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl<C: FpConfig<N>, const N: usize> MulAssign<&mut Self> for Fp<C, N> {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: &mut Self) {
        *self *= &*rhs;
    }
}

impl<C: FpConfig<N>, const N: usize> Div<&Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn div(mut self, rhs: &Self) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<C: FpConfig<N>, const N: usize> Div<Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn div(self, rhs: Self) -> Self::Output {
        self / &rhs
    }
}

impl<C: FpConfig<N>, const N: usize> Div<&mut Self> for Fp<C, N> {
    type Output = Self;

    #[inline(always)]
    fn div(self, rhs: &mut Self) -> Self::Output {
        self / &*rhs
    }
}

impl<C: FpConfig<N>, const N: usize> DivAssign<&Self> for Fp<C, N> {
    #[inline(always)]
    fn div_assign(&mut self, rhs: &Self) {
        *self *= rhs.inverse().expect("division by zero");
    }
}

impl<C: FpConfig<N>, const N: usize> DivAssign<Self> for Fp<C, N> {
    #[inline(always)]
    fn div_assign(&mut self, rhs: Self) {
        *self /= &rhs;
    }
}

impl<C: FpConfig<N>, const N: usize> DivAssign<&mut Self> for Fp<C, N> {
    #[inline(always)]
    fn div_assign(&mut self, rhs: &mut Self) {
        *self /= &*rhs;
    }
}

impl<C: FpConfig<N>, const N: usize> Sum<Self> for Fp<C, N> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), Add::add)
    }
}

impl<'a, C: FpConfig<N>, const N: usize> Sum<&'a Self> for Fp<C, N> {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), Add::add)
    }
}

impl<C: FpConfig<N>, const N: usize> Product<Self> for Fp<C, N> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Mul::mul)
    }
}

impl<'a, C: FpConfig<N>, const N: usize> Product<&'a Self> for Fp<C, N> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Mul::mul)
    }
}

#[inline(always)]
fn cmp_limbs<const N: usize>(a: &[u64; N], b: &[u64; N]) -> Ordering {
    if N == 4 {
        if a[3] != b[3] {
            return a[3].cmp(&b[3]);
        }
        if a[2] != b[2] {
            return a[2].cmp(&b[2]);
        }
        if a[1] != b[1] {
            return a[1].cmp(&b[1]);
        }
        return a[0].cmp(&b[0]);
    }

    for i in (0..N).rev() {
        match a[i].cmp(&b[i]) {
            Ordering::Equal => {}
            non_eq => return non_eq,
        }
    }
    Ordering::Equal
}

#[inline(always)]
fn add_raw<const N: usize>(a: [u64; N], b: [u64; N]) -> ([u64; N], bool) {
    if N == 4 {
        let (result, carry) = add_raw_4([a[0], a[1], a[2], a[3]], [b[0], b[1], b[2], b[3]]);
        let mut out = [0u64; N];
        out[..4].copy_from_slice(&result);
        return (out, carry);
    }

    let mut out = [0u64; N];
    let mut carry = false;

    for i in 0..N {
        (out[i], carry) = a[i].carrying_add(b[i], carry);
    }

    (out, carry)
}

#[inline(always)]
fn sub_raw<const N: usize>(a: [u64; N], b: [u64; N]) -> ([u64; N], bool) {
    if N == 4 {
        let (result, borrow) = sub_raw_4([a[0], a[1], a[2], a[3]], [b[0], b[1], b[2], b[3]]);
        let mut out = [0u64; N];
        out[..4].copy_from_slice(&result);
        return (out, borrow);
    }

    let mut out = [0u64; N];
    let mut borrow = false;

    for i in 0..N {
        (out[i], borrow) = a[i].borrowing_sub(b[i], borrow);
    }

    (out, borrow)
}

#[inline(always)]
fn add_raw_4(a: [u64; 4], b: [u64; 4]) -> ([u64; 4], bool) {
    let sum0 = a[0] as u128 + b[0] as u128;
    let z0 = sum0 as u64;
    let sum1 = a[1] as u128 + b[1] as u128 + (sum0 >> 64);
    let z1 = sum1 as u64;
    let sum2 = a[2] as u128 + b[2] as u128 + (sum1 >> 64);
    let z2 = sum2 as u64;
    let sum3 = a[3] as u128 + b[3] as u128 + (sum2 >> 64);
    let z3 = sum3 as u64;
    ([z0, z1, z2, z3], (sum3 >> 64) != 0)
}

#[inline(always)]
fn sub_raw_4(a: [u64; 4], b: [u64; 4]) -> ([u64; 4], bool) {
    let (z0, b0) = a[0].overflowing_sub(b[0]);
    let (z1, b1) = a[1].borrowing_sub(b[1], b0);
    let (z2, b2) = a[2].borrowing_sub(b[2], b1);
    let (z3, b3) = a[3].borrowing_sub(b[3], b2);
    ([z0, z1, z2, z3], b3)
}

#[derive(Clone, Copy)]
struct MontScratch<const N: usize> {
    limbs: [u64; N],
    extra: [u64; 2],
}

impl<const N: usize> MontScratch<N> {
    #[inline(always)]
    const fn zeroed() -> Self {
        Self {
            limbs: [0u64; N],
            extra: [0u64; 2],
        }
    }
}

impl<const N: usize> Index<usize> for MontScratch<N> {
    type Output = u64;

    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        if index < N {
            &self.limbs[index]
        } else {
            &self.extra[index - N]
        }
    }
}

impl<const N: usize> IndexMut<usize> for MontScratch<N> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index < N {
            &mut self.limbs[index]
        } else {
            &mut self.extra[index - N]
        }
    }
}

#[inline(always)]
fn mont_mul<const N: usize>(
    a: [u64; N],
    b: [u64; N],
    modulus: [u64; N],
    mont_inv: u64,
) -> [u64; N] {
    if N == 4 {
        let a4 = [a[0], a[1], a[2], a[3]];
        let b4 = [b[0], b[1], b[2], b[3]];
        let modulus4 = [modulus[0], modulus[1], modulus[2], modulus[3]];
        let result4 = mont_mul_4(a4, b4, modulus4, mont_inv);
        let mut out = [0u64; N];
        out[..4].copy_from_slice(&result4);
        return out;
    }

    mont_mul_generic(a, b, modulus, mont_inv)
}

#[inline(always)]
fn mont_reduce<const N: usize>(a: [u64; N], modulus: [u64; N], mont_inv: u64) -> [u64; N] {
    if N == 4 {
        let a4 = [a[0], a[1], a[2], a[3]];
        let modulus4 = [modulus[0], modulus[1], modulus[2], modulus[3]];
        let result4 = mont_reduce_4(a4, modulus4, mont_inv);
        let mut out = [0u64; N];
        out[..4].copy_from_slice(&result4);
        return out;
    }

    mont_reduce_generic(a, modulus, mont_inv)
}

#[inline(always)]
fn mont_mul_limb<const N: usize>(
    a: u64,
    b: [u64; N],
    modulus: [u64; N],
    mont_inv: u64,
) -> [u64; N] {
    if a == 0 {
        return [0u64; N];
    }

    let mut t = MontScratch::<N>::zeroed();
    let mut carry = 0u64;
    for j in 0..N {
        (t[j], carry) = a.carrying_mul_add(b[j], carry, t[j]);
    }
    add_to_high_word(&mut t, carry);
    mont_reduce_scratch(t, modulus, mont_inv)
}

#[inline(always)]
fn add_to_high_word_4(t: &mut [u64; 6], carry: u64) {
    let (sum, overflow) = t[4].carrying_add(carry, false);
    t[4] = sum;

    if overflow {
        let (sum, overflow) = t[5].carrying_add(1, false);
        t[5] = sum;
        debug_assert!(!overflow);
    }
}

#[inline(always)]
fn mont_mul_4(a: [u64; 4], b: [u64; 4], modulus: [u64; 4], mont_inv: u64) -> [u64; 4] {
    let mut t = [0u64; 6];

    for i in 0..4 {
        let mut carry = 0u64;
        for j in 0..4 {
            (t[j], carry) = a[i].carrying_mul_add(b[j], carry, t[j]);
        }
        add_to_high_word_4(&mut t, carry);

        let m = t[0].wrapping_mul(mont_inv);
        carry = 0;
        for j in 0..4 {
            (t[j], carry) = m.carrying_mul_add(modulus[j], carry, t[j]);
        }
        debug_assert_eq!(t[0], 0);
        add_to_high_word_4(&mut t, carry);

        t[0] = t[1];
        t[1] = t[2];
        t[2] = t[3];
        t[3] = t[4];
        t[4] = t[5];
        t[5] = 0;
    }

    let mut result = [t[0], t[1], t[2], t[3]];
    debug_assert_eq!(t[4], 0);

    if cmp_limbs(&result, &modulus) != Ordering::Less {
        result = sub_raw(result, modulus).0;
    }

    result
}

#[inline(always)]
fn mont_reduce_4(a: [u64; 4], modulus: [u64; 4], mont_inv: u64) -> [u64; 4] {
    let mut t = [a[0], a[1], a[2], a[3], 0, 0];

    for _ in 0..4 {
        let m = t[0].wrapping_mul(mont_inv);
        let mut carry = 0u64;
        for j in 0..4 {
            (t[j], carry) = m.carrying_mul_add(modulus[j], carry, t[j]);
        }
        debug_assert_eq!(t[0], 0);
        add_to_high_word_4(&mut t, carry);

        t[0] = t[1];
        t[1] = t[2];
        t[2] = t[3];
        t[3] = t[4];
        t[4] = t[5];
        t[5] = 0;
    }

    let mut result = [t[0], t[1], t[2], t[3]];
    debug_assert_eq!(t[4], 0);

    if cmp_limbs(&result, &modulus) != Ordering::Less {
        result = sub_raw(result, modulus).0;
    }

    result
}

#[inline(always)]
fn add_to_high_word<const N: usize>(t: &mut MontScratch<N>, mut carry: u64) {
    let mut k = N;
    while carry != 0 {
        let (sum, overflow) = t[k].carrying_add(carry, false);
        t[k] = sum;
        carry = overflow as u64;
        k += 1;
    }
}

#[inline(always)]
fn shift_scratch<const N: usize>(t: &mut MontScratch<N>) {
    for j in 0..=N {
        t[j] = t[j + 1];
    }
    t[N + 1] = 0;
}

#[inline(always)]
fn mont_reduce_scratch<const N: usize>(
    mut t: MontScratch<N>,
    modulus: [u64; N],
    mont_inv: u64,
) -> [u64; N] {
    for _ in 0..N {
        let m = t[0].wrapping_mul(mont_inv);
        let mut carry = 0u64;
        for j in 0..N {
            (t[j], carry) = m.carrying_mul_add(modulus[j], carry, t[j]);
        }
        debug_assert_eq!(t[0], 0);
        add_to_high_word(&mut t, carry);
        shift_scratch(&mut t);
    }

    let mut result = [0u64; N];
    result.copy_from_slice(&t.limbs);
    debug_assert_eq!(t[N], 0);

    if cmp_limbs(&result, &modulus) != Ordering::Less {
        result = sub_raw(result, modulus).0;
    }

    result
}

#[inline(always)]
fn mont_reduce_generic<const N: usize>(a: [u64; N], modulus: [u64; N], mont_inv: u64) -> [u64; N] {
    let mut t = MontScratch::<N>::zeroed();
    t.limbs = a;
    mont_reduce_scratch(t, modulus, mont_inv)
}

#[inline(always)]
fn mont_mul_generic<const N: usize>(
    a: [u64; N],
    b: [u64; N],
    modulus: [u64; N],
    mont_inv: u64,
) -> [u64; N] {
    let mut t = MontScratch::<N>::zeroed();

    for i in 0..N {
        let mut carry = 0u64;
        for j in 0..N {
            (t[j], carry) = a[i].carrying_mul_add(b[j], carry, t[j]);
        }
        add_to_high_word(&mut t, carry);

        let m = t[0].wrapping_mul(mont_inv);
        carry = 0;
        for j in 0..N {
            (t[j], carry) = m.carrying_mul_add(modulus[j], carry, t[j]);
        }
        debug_assert_eq!(t[0], 0);
        add_to_high_word(&mut t, carry);

        shift_scratch(&mut t);
    }

    let mut result = [0u64; N];
    result.copy_from_slice(&t.limbs);
    debug_assert_eq!(t[N], 0);

    if cmp_limbs(&result, &modulus) != Ordering::Less {
        result = sub_raw(result, modulus).0;
    }

    result
}

fn reduce_limbs<const N: usize>(x: [u64; N], modulus: [u64; N]) -> [u64; N] {
    if cmp_limbs(&x, &modulus) == Ordering::Less {
        return x;
    }

    let mut reduced = x;
    for _ in 0..fast_reduce_subtractions(&modulus) {
        let (next, borrow) = sub_raw(reduced, modulus);
        if borrow {
            break;
        }
        if cmp_limbs(&next, &modulus) == Ordering::Less {
            return next;
        }
        reduced = next;
    }

    const_reduce(x, modulus)
}

#[inline(always)]
fn fast_reduce_subtractions<const N: usize>(modulus: &[u64; N]) -> usize {
    let modulus_bits = limb_bit_len(modulus);
    let unused_bits = N * 64 - modulus_bits;
    let mut attempts = 1usize;
    let mut i = 0;
    while i <= unused_bits && attempts < 64 {
        attempts <<= 1;
        i += 1;
    }
    attempts
}

#[inline(always)]
fn limb_bit_len<const N: usize>(limbs: &[u64; N]) -> usize {
    let mut i = N;
    while i > 0 {
        i -= 1;
        let limb = limbs[i];
        if limb != 0 {
            return i * 64 + (64 - limb.leading_zeros() as usize);
        }
    }
    0
}

#[inline(always)]
fn high_limbs_are_zero<const N: usize>(limbs: &[u64; N]) -> bool {
    let mut i = 1;
    while i < N {
        if limbs[i] != 0 {
            return false;
        }
        i += 1;
    }
    true
}

fn mask_unused_limb_bits<const N: usize>(limbs: &mut [u64; N], modulus_bits: usize) {
    let unused_bits = 64 * N - modulus_bits;
    if unused_bits > 0 {
        limbs[N - 1] &= u64::MAX >> unused_bits;
    }
}

const fn one_limbs<const N: usize>() -> [u64; N] {
    let mut limbs = [0u64; N];
    limbs[0] = 1;
    limbs
}

const fn mont_inv<const N: usize>(modulus: [u64; N]) -> u64 {
    let p0 = modulus[0];
    let mut inv = 1u64;
    let mut i = 0;
    while i < 6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(p0.wrapping_mul(inv)));
        i += 1;
    }
    inv.wrapping_neg()
}

const fn const_cmp<const N: usize>(a: [u64; N], b: [u64; N]) -> i8 {
    let mut i = N;
    while i > 0 {
        i -= 1;
        if a[i] < b[i] {
            return -1;
        }
        if a[i] > b[i] {
            return 1;
        }
    }
    0
}

const fn const_add_raw<const N: usize>(a: [u64; N], b: [u64; N]) -> ([u64; N], bool) {
    let mut out = [0u64; N];
    let mut carry = 0u128;
    let mut i = 0;
    while i < N {
        let sum = a[i] as u128 + b[i] as u128 + carry;
        out[i] = sum as u64;
        carry = sum >> 64;
        i += 1;
    }
    (out, carry != 0)
}

const fn const_sub_raw<const N: usize>(a: [u64; N], b: [u64; N]) -> ([u64; N], bool) {
    let mut out = [0u64; N];
    let mut borrow = 0u128;
    let mut i = 0;
    while i < N {
        let subtrahend = b[i] as u128 + borrow;
        let ai = a[i] as u128;
        if ai >= subtrahend {
            out[i] = (ai - subtrahend) as u64;
            borrow = 0;
        } else {
            out[i] = ((1u128 << 64) + ai - subtrahend) as u64;
            borrow = 1;
        }
        i += 1;
    }
    (out, borrow != 0)
}

const fn const_add_mod<const N: usize>(a: [u64; N], b: [u64; N], modulus: [u64; N]) -> [u64; N] {
    let (mut sum, carry) = const_add_raw(a, b);
    if carry || const_cmp(sum, modulus) >= 0 {
        sum = const_sub_raw(sum, modulus).0;
    }
    sum
}

const fn const_double_mod<const N: usize>(a: [u64; N], modulus: [u64; N]) -> [u64; N] {
    const_add_mod(a, a, modulus)
}

const fn const_get_bit<const N: usize>(a: [u64; N], bit: usize) -> bool {
    let limb = bit / 64;
    if limb >= N {
        false
    } else {
        (a[limb] & (1u64 << (bit % 64))) != 0
    }
}

const fn const_reduce<const N: usize>(x: [u64; N], modulus: [u64; N]) -> [u64; N] {
    let mut result = [0u64; N];
    let mut i = N * 64;
    while i > 0 {
        i -= 1;
        result = const_double_mod(result, modulus);
        if const_get_bit(x, i) {
            result = const_add_mod(result, one_limbs(), modulus);
        }
    }
    result
}

const fn const_mul_mod<const N: usize>(a: [u64; N], b: [u64; N], modulus: [u64; N]) -> [u64; N] {
    let mut result = [0u64; N];
    let mut base = const_reduce(a, modulus);
    let mut i = 0;
    while i < N * 64 {
        if const_get_bit(b, i) {
            result = const_add_mod(result, base, modulus);
        }
        base = const_double_mod(base, modulus);
        i += 1;
    }
    result
}

const fn const_pow_mod<const N: usize>(
    base: [u64; N],
    exp: [u64; N],
    modulus: [u64; N],
) -> [u64; N] {
    let mut result = one_limbs();
    let mut power = const_reduce(base, modulus);
    let mut i = 0;
    while i < N * 64 {
        if const_get_bit(exp, i) {
            result = const_mul_mod(result, power, modulus);
        }
        power = const_mul_mod(power, power, modulus);
        i += 1;
    }
    result
}

const fn mont_r<const N: usize>(modulus: [u64; N]) -> [u64; N] {
    let mut result = one_limbs();
    let mut i = 0;
    while i < N * 64 {
        result = const_double_mod(result, modulus);
        i += 1;
    }
    result
}

const fn mont_r2<const N: usize>(modulus: [u64; N]) -> [u64; N] {
    let mut result = one_limbs();
    let mut i = 0;
    while i < N * 128 {
        result = const_double_mod(result, modulus);
        i += 1;
    }
    result
}

const fn to_montgomery<const N: usize>(value: [u64; N], modulus: [u64; N]) -> [u64; N] {
    const_mul_mod(value, mont_r(modulus), modulus)
}

const fn mont_neg_one<const N: usize>(modulus: [u64; N]) -> [u64; N] {
    let one = mont_r(modulus);
    if const_cmp(one, [0u64; N]) == 0 {
        one
    } else {
        const_sub_raw(modulus, one).0
    }
}

const fn modulus_minus_two<const N: usize>(modulus: [u64; N]) -> [u64; N] {
    let mut two = [0u64; N];
    two[0] = 2;
    const_sub_raw(modulus, two).0
}

const fn div2<const N: usize>(mut value: [u64; N]) -> [u64; N] {
    let mut carry = 0u64;
    let mut i = N;
    while i > 0 {
        i -= 1;
        let next_carry = value[i] << 63;
        value[i] = (value[i] >> 1) | carry;
        carry = next_carry;
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq as ArkBls12381Fq;
    use ark_ff::{BigInt as ArkBigInt, MontBackend, MontConfig, PrimeField};
    use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

    #[derive(MontConfig)]
    #[modulus = "57896044618658097711785492504343953926634992332820282019728792003956564819949"]
    #[generator = "2"]
    struct ArkF25519Config;

    type ArkFp = ark_ff::Fp256<MontBackend<ArkF25519Config, 4>>;

    #[derive(Debug)]
    struct TestBls12381FqConfig;

    impl FpConfig<6> for TestBls12381FqConfig {
        const MODULUS: [u64; 6] = <ArkBls12381Fq as PrimeField>::MODULUS.0;
        const GENERATOR: [u64; 6] = [2, 0, 0, 0, 0, 0];
    }

    type TestBls12381Fq = Fp<TestBls12381FqConfig, 6>;

    fn ark_from_limbs(limbs: [u64; NUM_LIMBS]) -> ArkFp {
        ArkFp::from_le_bytes_mod_order(&limbs_to_fixed_le_bytes(limbs))
    }

    fn limbs_to_fixed_le_bytes(limbs: [u64; NUM_LIMBS]) -> [u8; 32] {
        let mut bytes = [0u8; 32];
        for (i, limb) in limbs.iter().enumerate() {
            bytes[i * 8..(i + 1) * 8].copy_from_slice(&limb.to_le_bytes());
        }
        bytes
    }

    fn bls_limbs_to_fixed_le_bytes(limbs: [u64; 6]) -> [u8; 48] {
        let mut bytes = [0u8; 48];
        for (i, limb) in limbs.iter().enumerate() {
            bytes[i * 8..(i + 1) * 8].copy_from_slice(&limb.to_le_bytes());
        }
        bytes
    }

    fn ark_bls_from_limbs(limbs: [u64; 6]) -> ArkBls12381Fq {
        ArkBls12381Fq::from_le_bytes_mod_order(&bls_limbs_to_fixed_le_bytes(limbs))
    }

    fn ark_limbs(x: ArkFp) -> [u64; NUM_LIMBS] {
        x.into_bigint().0
    }

    fn next_limbs(state: &mut u64) -> [u64; NUM_LIMBS] {
        let mut limbs = [0u64; NUM_LIMBS];
        for limb in &mut limbs {
            *state = state
                .wrapping_mul(0xda94_2042_e4dd_58b5)
                .wrapping_add(0x9e37_79b9_7f4a_7c15);
            *limb = *state;
        }
        limbs
    }

    fn assert_prime_field<F: ark_ff::PrimeField>() {}

    #[test]
    fn implements_ark_prime_field_traits() {
        assert_prime_field::<Fp25519>();
    }

    #[test]
    fn constants_are_generated_at_const_layer() {
        assert_eq!(Fp25519::MONT_INV, 0x86bca1af286bca1b);
        assert_eq!(Fp25519::R2, [1444, 0, 0, 0]);
        assert_eq!(Fp25519::MODULUS, MODULUS_LIMBS);
        assert_eq!(Fp25519::ONE.to_limbs(), [1, 0, 0, 0]);
        assert_eq!(
            Fp25519::TWO_ADIC_ROOT_OF_UNITY.to_limbs(),
            ArkFp::TWO_ADIC_ROOT_OF_UNITY.into_bigint().0
        );
    }

    #[test]
    fn constructors_round_trip_and_reduce() {
        assert_eq!(Fp25519::zero().to_limbs(), [0, 0, 0, 0]);
        assert_eq!(Fp25519::one().to_limbs(), [1, 0, 0, 0]);
        assert_eq!(Fp25519::from_u64(42).to_limbs(), [42, 0, 0, 0]);
        assert_eq!(Fp25519::new(MODULUS_LIMBS).to_limbs(), [0, 0, 0, 0]);

        let mut p_plus_five = MODULUS_LIMBS;
        p_plus_five[0] += 5;
        assert_eq!(Fp25519::new(p_plus_five).to_limbs(), [5, 0, 0, 0]);
    }

    #[test]
    fn n6_constructors_match_ark_ff() {
        assert_eq!(fast_reduce_subtractions(&TestBls12381FqConfig::MODULUS), 16);

        let vectors = [
            [0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [u64::MAX, 0, 0, 0, 0, 0],
            [u64::MAX, u64::MAX, u64::MAX, u64::MAX, u64::MAX, u64::MAX],
            [
                0x0123456789abcdef,
                0xfedcba9876543210,
                0x0f0f0f0f0f0f0f0f,
                0xf0f0f0f0f0f0f0f0,
                0xaaaaaaaa55555555,
                0xffffffffffffffff,
            ],
        ];

        for limbs in vectors {
            assert_eq!(
                TestBls12381Fq::new(limbs).to_limbs(),
                ark_bls_from_limbs(limbs).into_bigint().0
            );
        }

        for value in [0, 1, 42, u64::MAX] {
            assert_eq!(
                TestBls12381Fq::from_u64(value).to_limbs(),
                ArkBls12381Fq::from(value).into_bigint().0
            );
        }
    }

    #[test]
    fn small_arithmetic_examples() {
        let two = Fp25519::from_u64(2);
        let three = Fp25519::from_u64(3);

        assert_eq!((two + three).to_limbs(), [5, 0, 0, 0]);
        assert_eq!((two * three).to_limbs(), [6, 0, 0, 0]);
        assert_eq!((three - two).to_limbs(), [1, 0, 0, 0]);

        assert_eq!(
            (two - three).to_limbs(),
            [
                0xffffffffffffffec,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0x7fffffffffffffff,
            ],
        );
    }

    #[test]
    fn arithmetic_matches_ark_ff_on_vectors() {
        let vectors = [
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [2, 0, 0, 0],
            [19, 0, 0, 0],
            [u64::MAX, 0, 0, 0],
            [0, 1, 0, 0],
            [
                0x0123456789abcdef,
                0xfedcba9876543210,
                0x0f0f0f0f0f0f0f0f,
                0x7f0f0f0f0f0f0f0f,
            ],
            [
                0xfffffffffffffff0,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0x7fffffffffffffff,
            ],
            [
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
            ],
        ];

        for a_limbs in vectors {
            for b_limbs in vectors {
                let a = Fp25519::new(a_limbs);
                let b = Fp25519::new(b_limbs);
                let ark_a = ark_from_limbs(a_limbs);
                let ark_b = ark_from_limbs(b_limbs);

                assert_eq!((a + b).to_limbs(), ark_limbs(ark_a + ark_b));
                assert_eq!((a - b).to_limbs(), ark_limbs(ark_a - ark_b));
                assert_eq!((a * b).to_limbs(), ark_limbs(ark_a * ark_b));
            }
        }
    }

    #[test]
    fn multiplication_matches_ark_ff_on_deterministic_inputs() {
        let mut state = 0x1319_8a2e_0370_7344u64;

        for _ in 0..256 {
            let a_limbs = next_limbs(&mut state);
            let b_limbs = next_limbs(&mut state);
            let a = Fp25519::new(a_limbs);
            let b = Fp25519::new(b_limbs);
            let ark_a = ark_from_limbs(a_limbs);
            let ark_b = ark_from_limbs(b_limbs);

            assert_eq!((a * b).to_limbs(), ark_limbs(ark_a * ark_b));
        }
    }

    #[test]
    fn mul_batch_matches_scalar_mul() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let lhs = (0..64)
            .map(|_| Fp25519::new(next_limbs(&mut state)))
            .collect::<Vec<_>>();
        let rhs = (0..64)
            .map(|_| Fp25519::new(next_limbs(&mut state)))
            .collect::<Vec<_>>();
        let mut out = vec![Fp25519::zero(); lhs.len()];

        Fp25519::mul_batch(&lhs, &rhs, &mut out);

        for i in 0..lhs.len() {
            assert_eq!(out[i], lhs[i] * rhs[i]);
        }
    }

    #[test]
    #[should_panic]
    fn mul_batch_panics_on_length_mismatch() {
        let lhs = [Fp25519::one()];
        let rhs = [Fp25519::one(), Fp25519::one()];
        let mut out = [Fp25519::zero()];

        Fp25519::mul_batch(&lhs, &rhs, &mut out);
    }

    #[test]
    fn canonical_limbs_can_be_loaded_by_ark_ff() {
        let x = Fp25519::new([
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
        ]);

        assert_eq!(
            ArkFp::from_bigint(ArkBigInt(x.to_limbs())).unwrap(),
            ark_from_limbs([
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
            ]),
        );
    }

    #[test]
    fn bigint_and_serialization_round_trip() {
        let x = Fp25519::new([
            0x0123456789abcdef,
            0xfedcba9876543210,
            0x0f0f0f0f0f0f0f0f,
            0x7f0f0f0f0f0f0f0f,
        ]);

        assert_eq!(Fp25519::from_bigint(x.into_bigint()), Some(x));

        let mut bytes = Vec::new();
        x.serialize_compressed(&mut bytes).unwrap();
        let y = Fp25519::deserialize_compressed(bytes.as_slice()).unwrap();
        assert_eq!(x, y);
    }
}
