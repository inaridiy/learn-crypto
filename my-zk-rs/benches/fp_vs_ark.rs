use std::{hint::black_box, ops::Add, ops::Mul, ops::Sub, time::Duration};

use ark_bls12_381::{Fq as ArkBls12381Fq, Fr as ArkBls12381Fr};
use ark_bn254::{Fq as ArkBn254Fq, Fr as ArkBn254Fr};
use ark_ff::{BigInt as ArkBigInt, MontBackend, MontConfig, PrimeField};
use ark_pallas::Fq as ArkPallasFq;
use criterion::{Criterion, criterion_group, criterion_main};
use my_zk_rs::{F25519Config, Fp as MyFp, FpConfig};

#[derive(MontConfig)]
#[modulus = "57896044618658097711785492504343953926634992332820282019728792003956564819949"]
#[generator = "2"]
struct ArkF25519Config;

type ArkFp25519 = ark_ff::Fp256<MontBackend<ArkF25519Config, 4>>;

#[derive(Debug)]
struct Bn254FrConfig;

impl FpConfig<4> for Bn254FrConfig {
    const MODULUS: [u64; 4] = <ArkBn254Fr as PrimeField>::MODULUS.0;
    const GENERATOR: [u64; 4] = [5, 0, 0, 0];
}

#[derive(Debug)]
struct Bn254FqConfig;

impl FpConfig<4> for Bn254FqConfig {
    const MODULUS: [u64; 4] = <ArkBn254Fq as PrimeField>::MODULUS.0;
    const GENERATOR: [u64; 4] = [5, 0, 0, 0];
}

#[derive(Debug)]
struct Bls12381FrConfig;

impl FpConfig<4> for Bls12381FrConfig {
    const MODULUS: [u64; 4] = <ArkBls12381Fr as PrimeField>::MODULUS.0;
    const GENERATOR: [u64; 4] = [7, 0, 0, 0];
}

#[derive(Debug)]
struct Bls12381FqConfig;

impl FpConfig<6> for Bls12381FqConfig {
    const MODULUS: [u64; 6] = <ArkBls12381Fq as PrimeField>::MODULUS.0;
    const GENERATOR: [u64; 6] = [2, 0, 0, 0, 0, 0];
}

#[derive(Debug)]
struct PallasFqConfig;

impl FpConfig<4> for PallasFqConfig {
    const MODULUS: [u64; 4] = <ArkPallasFq as PrimeField>::MODULUS.0;
    const GENERATOR: [u64; 4] = [5, 0, 0, 0];
}

const ARK_LABEL: &str = if cfg!(feature = "ark-asm") {
    "ark_ff_asm"
} else {
    "ark_ff"
};

const SAMPLES: usize = 1024;

fn sample_limbs<C: FpConfig<N>, const N: usize>() -> Vec<[u64; N]> {
    let mut state = 0x243f_6a88_85a3_08d3u64 ^ ((N as u64) << 32);
    let mut samples = Vec::with_capacity(SAMPLES);

    for _ in 0..SAMPLES {
        let mut limbs = [0u64; N];
        for limb in &mut limbs {
            state = state
                .wrapping_mul(0xda94_2042_e4dd_58b5)
                .wrapping_add(0x9e37_79b9_7f4a_7c15);
            *limb = state;
        }

        let reduced = MyFp::<C, N>::new(limbs).to_limbs();
        samples.push(if reduced == [0u64; N] {
            let mut one = [0u64; N];
            one[0] = 1;
            one
        } else {
            reduced
        });
    }

    samples
}

fn make_input_pairs<C, A, const N: usize>() -> (Vec<MyFp<C, N>>, Vec<MyFp<C, N>>, Vec<A>, Vec<A>)
where
    C: FpConfig<N>,
    A: PrimeField<BigInt = ArkBigInt<N>> + Copy,
{
    let limbs = sample_limbs::<C, N>();
    let my_lhs = limbs.iter().copied().map(MyFp::<C, N>::new).collect();
    let ark_lhs = limbs
        .iter()
        .copied()
        .map(|limbs| A::from_bigint(ArkBigInt(limbs)).unwrap())
        .collect::<Vec<_>>();
    let my_rhs = limbs.iter().copied().rev().map(MyFp::<C, N>::new).collect();
    let ark_rhs = limbs
        .into_iter()
        .rev()
        .map(|limbs| A::from_bigint(ArkBigInt(limbs)).unwrap())
        .collect();

    (my_lhs, my_rhs, ark_lhs, ark_rhs)
}

fn assert_matches_ark<C, A, const N: usize>(
    my_lhs: &[MyFp<C, N>],
    my_rhs: &[MyFp<C, N>],
    ark_lhs: &[A],
    ark_rhs: &[A],
) where
    C: FpConfig<N>,
    A: PrimeField<BigInt = ArkBigInt<N>>
        + Copy
        + Add<Output = A>
        + Sub<Output = A>
        + Mul<Output = A>,
{
    for i in 0..32 {
        assert_eq!(
            (my_lhs[i] + my_rhs[i]).to_limbs(),
            (ark_lhs[i] + ark_rhs[i]).into_bigint().0
        );
        assert_eq!(
            (my_lhs[i] - my_rhs[i]).to_limbs(),
            (ark_lhs[i] - ark_rhs[i]).into_bigint().0
        );
        assert_eq!(
            (my_lhs[i] * my_rhs[i]).to_limbs(),
            (ark_lhs[i] * ark_rhs[i]).into_bigint().0
        );
    }
}

fn bench_field<C, A, const N: usize>(c: &mut Criterion, name: &str)
where
    C: FpConfig<N>,
    A: PrimeField<BigInt = ArkBigInt<N>>
        + Copy
        + Add<Output = A>
        + Sub<Output = A>
        + Mul<Output = A>,
{
    let (my_lhs, my_rhs, ark_lhs, ark_rhs) = make_input_pairs::<C, A, N>();
    assert_matches_ark(&my_lhs, &my_rhs, &ark_lhs, &ark_rhs);

    let mut group = c.benchmark_group(format!("add/{name}"));
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_fp", |b| {
        b.iter(|| {
            let mut acc = MyFp::<C, N>::zero();
            for x in &my_lhs {
                acc = black_box(acc) + black_box(*x);
            }
            black_box(acc)
        })
    });

    group.bench_function(ARK_LABEL, |b| {
        b.iter(|| {
            let mut acc = A::zero();
            for x in &ark_lhs {
                acc = black_box(acc) + black_box(*x);
            }
            black_box(acc)
        })
    });
    group.finish();

    let mut group = c.benchmark_group(format!("sub/{name}"));
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_fp", |b| {
        b.iter(|| {
            let mut acc = MyFp::<C, N>::one();
            for x in &my_lhs {
                acc = black_box(acc) - black_box(*x);
            }
            black_box(acc)
        })
    });

    group.bench_function(ARK_LABEL, |b| {
        b.iter(|| {
            let mut acc = A::one();
            for x in &ark_lhs {
                acc = black_box(acc) - black_box(*x);
            }
            black_box(acc)
        })
    });
    group.finish();

    let mut group = c.benchmark_group(format!("mul/{name}"));
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_fp", |b| {
        b.iter(|| {
            let mut acc = MyFp::<C, N>::one();
            for x in &my_lhs {
                acc = black_box(acc) * black_box(*x);
            }
            black_box(acc)
        })
    });

    group.bench_function(ARK_LABEL, |b| {
        b.iter(|| {
            let mut acc = A::one();
            for x in &ark_lhs {
                acc = black_box(acc) * black_box(*x);
            }
            black_box(acc)
        })
    });
    group.finish();

    let mut group = c.benchmark_group(format!("batch_mul/{name}"));
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(1));
    group.sample_size(30);

    group.bench_function("my_scalar", |b| {
        let mut out = vec![MyFp::<C, N>::zero(); my_lhs.len()];
        b.iter(|| {
            MyFp::<C, N>::mul_batch(black_box(&my_lhs), black_box(&my_rhs), black_box(&mut out));
            black_box(out[0])
        })
    });

    group.bench_function(
        if cfg!(feature = "ark-asm") {
            "ark_ff_asm_loop"
        } else {
            "ark_ff_loop"
        },
        |b| {
            let mut out = vec![A::zero(); ark_lhs.len()];
            b.iter(|| {
                for ((out, lhs), rhs) in out.iter_mut().zip(&ark_lhs).zip(&ark_rhs) {
                    *out = black_box(*lhs) * black_box(*rhs);
                }
                black_box(out[0])
            })
        },
    );
    group.finish();
}

fn bench_common_moduli(c: &mut Criterion) {
    bench_field::<F25519Config, ArkFp25519, 4>(c, "fp25519");
    bench_field::<Bn254FrConfig, ArkBn254Fr, 4>(c, "bn254_fr");
    bench_field::<Bn254FqConfig, ArkBn254Fq, 4>(c, "bn254_fq");
    bench_field::<Bls12381FrConfig, ArkBls12381Fr, 4>(c, "bls12_381_fr");
    bench_field::<Bls12381FqConfig, ArkBls12381Fq, 6>(c, "bls12_381_fq");
    bench_field::<PallasFqConfig, ArkPallasFq, 4>(c, "pallas_fq");
}

criterion_group!(benches, bench_common_moduli);
criterion_main!(benches);
