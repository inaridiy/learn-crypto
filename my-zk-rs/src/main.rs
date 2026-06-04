#![feature(generic_const_exprs)]
#![feature(min_adt_const_params)]

use std::{cmp::Ordering, ops::Index};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct BigInt<const N: usize>([u64; N]);

impl<const N: usize> Ord for BigInt<N> {
    fn cmp(&self, rhs: &Self) -> Ordering {
        for i in (0..N).rev() {
            match self[i].cmp(&rhs[i]) {
                Ordering::Equal => {}
                non_eq => return non_eq,
            }
        }

        Ordering::Equal
    }
}

impl<const N: usize> PartialOrd for BigInt<N> {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}
impl<const N: usize> From<BigInt<N>> for [u64; N] {
    fn from(value: BigInt<N>) -> Self {
        value.0
    }
}

impl<const N: usize> From<[u64; N]> for BigInt<N> {
    fn from(value: [u64; N]) -> Self {
        BigInt(value)
    }
}

impl<const N: usize> Index<usize> for BigInt<N> {
    type Output = u64;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

#[inline(always)]
fn add_raw<const N: usize>(a: BigInt<N>, b: BigInt<N>) -> (BigInt<N>, bool) {
    let mut out = [0u64; N];
    let mut carry = false;

    for i in 0..N {
        (out[i], carry) = a[i].carrying_add(b[i], carry);
    }

    (BigInt(out), carry)
}

#[inline(always)]
fn sub_raw<const N: usize>(a: BigInt<N>, b: BigInt<N>) -> (BigInt<N>, bool) {
    let mut out = [0u64; N];
    let mut borrow = false;

    for i in 0..N {
        (out[i], borrow) = a[i].borrowing_sub(b[i], borrow);
    }

    (BigInt(out), borrow)
}

#[inline(always)]
fn mul_raw_full<const N: usize, const O: usize>(a: BigInt<N>, b: BigInt<N>) -> BigInt<O> {
    let mut out = [0u64; O];

    for i in 0..N {
        let mut carry = 0u64;

        for j in 0..N {
            let k = i + j;
            (out[k], carry) = a[i].carrying_mul_add(b[j], carry, out[k]);
        }

        let mut k = i + N;
        while carry != 0 {
            let (sum, overflow) = out[k].carrying_add(carry, false);
            out[k] = sum;
            carry = overflow as u64;
            k += 1;
        }
    }

    BigInt(out)
}

//------

const LIMBS: usize = 4;
const WIDE_LIMBS: usize = LIMBS * 2 + 1;

// P = 2^{255} - 19
const MODULUS: BigInt<LIMBS> = BigInt([
    0xffffffffffffffed,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0x7fffffffffffffff,
]);

// -P^{-1} mod 2^64
const MONT_INV: u64 = 0x86bca1af286bca1b;

// R2 = R^2 mod P
const R2: BigInt<LIMBS> = BigInt([1444, 0, 0, 0]);

fn mont_reduction(mut t: BigInt<WIDE_LIMBS>) -> BigInt<LIMBS> {
    for i in 0..LIMBS {
        // m = t[i] * (-P^{-1}) mod 2^{64}
        // この m を使うと、下で t に m * MODULUS * 2^{64*i} を足したとき
        // t[i] がちょうど 0 になる（= R で割れる桁が 1 つ増える）。
        let m = t[i].wrapping_mul(MONT_INV);

        // t += m * MODULUS << (64 * i)
        let mut carry: u64 = 0;
        for j in 0..LIMBS {
            let k = i + j;
            (t.0[k], carry) = m.carrying_mul_add(MODULUS[j], carry, t[k]);
        }

        // 内側ループから出た carry を上位ワードへ伝播させる。
        // 足し込みは t[i + LIMBS] 以降に対して行う（t[i] ではない）。
        let mut k = i + LIMBS;
        while carry != 0 {
            let (sum, overflow) = t[k].carrying_add(carry, false);
            t.0[k] = sum;
            carry = overflow as u64;

            k += 1;
        }
    }

    // 上位 LIMBS ワードが T * R^{-1} (mod P) の候補。値の範囲は [0, 2P)。
    let mut out = [0u64; LIMBS];
    out.copy_from_slice(&t.0[LIMBS..LIMBS * 2]);
    let mut result = BigInt(out);

    // 最終ワードの桁上がり。入力が正しい範囲なら 0 になるはず。
    let high = t[LIMBS * 2];
    debug_assert_eq!(high, 0);

    // 候補が [P, 2P) なら 1 回だけ MODULUS を引いて [0, P) に収める。
    if high != 0 || result >= MODULUS {
        let (reduced, _borrow) = sub_raw(result, MODULUS);
        result = reduced;
    }

    result
}

fn mont_mul(a: BigInt<LIMBS>, b: BigInt<LIMBS>) -> BigInt<LIMBS> {
    let t = mul_raw_full(a, b);
    mont_reduction(t)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp {
    limbs: BigInt<LIMBS>,
}

impl Fp {
    pub fn new(x: [u64; LIMBS]) -> Self {
        let mut x = BigInt(x);
        while x >= MODULUS {
            (x, _) = sub_raw(x, MODULUS);
        }
        Self {
            limbs: mont_mul(x, R2),
        }
    }

    pub fn add(self, rhs: Self) -> Self {
        let (mut sum, carry) = add_raw(self.limbs, rhs.limbs);

        if carry {
            sum = sub_raw(sum, MODULUS).0;
        }

        if sum >= MODULUS {
            sum = sub_raw(sum, MODULUS).0;
        }

        Self { limbs: sum }
    }

    /// Montgomery 表現同士の減算。
    pub fn sub(self, rhs: Self) -> Self {
        let (mut diff, borrow) = sub_raw(self.limbs, rhs.limbs);

        if borrow {
            diff = add_raw(diff, MODULUS).0;
        }

        Self { limbs: diff }
    }

    /// Montgomery 表現同士の乗算。
    pub fn mul(self, rhs: Self) -> Self {
        Self {
            limbs: mont_mul(self.limbs, rhs.limbs),
        }
    }
}

fn main() {
    let R = (u64::MAX as u128) + 1; // 2^{64}
}
