//! Zero-knowledge sumcheck。
//!
//! 平文の [`SumcheckProof`] はラウンド多項式 $s_i$ を
//! そのまま送るが、$s_i$ は witness の情報を漏らす。ここでは Spartan の参照実装に倣い、
//! prover は各ラウンドで
//!
//! 1. $s_i$ の係数ベクトルへの Pedersen commitment $C_{s_i}$、
//! 2. challenge $r_i$ での評価 $s_i(r_i)$ への commitment $C_{e_i}$(次ラウンドの claim)、
//! 3. 「$s_i(0) + s_i(1) = e_{i-1}$」と「$s_i(r_i) = e_i$」を random linear combination で
//!    一本にまとめた dot-product argument([`LinearInnerProductProof`])
//!
//! を送る。係数ベクトルと $(2, 1, \ldots, 1)$ の内積が $s_i(0) + s_i(1)$、
//! $(1, r_i, r_i^2, \ldots)$ との内積が $s_i(r_i)$ なので、どちらも内積関係になる。
//! verifier は claim を commitment としてしか持たず、最終 claim も commitment
//! $C_{e_{n-1}}$ として返す。上位のプロトコルは、それを別の sigma protocol で処理する。

use ark_ec::{
    hashing::{HashToCurve, HashToCurveError},
    CurveGroup,
};
use ark_ff::{Field, UniformRand};
use ark_std::rand::{CryptoRng, Rng};

use crate::{
    primitive::{MultilinearPoly, ScalarPedersen, Transcript, UniPoly, VectorPedersen},
    protocol::{sigma::LinearInnerProductProof, sumcheck::SumcheckProof},
};

/// ZK sumcheck の commitment key。
///
/// ラウンド多項式の係数(長さ `degree + 1`)を commit する vector key と、
/// claim を commit する scalar key。claim の commitment は上位のプロトコルと
/// 準同型に組み合わせるので、scalar key は呼び出し側から受け取り、blind generator も共有する。
#[derive(Clone, Debug)]
pub struct ZkSumcheckKey<G: CurveGroup> {
    pub coeffs: VectorPedersen<G>,
    pub scalar: ScalarPedersen<G>,
}

impl<G: CurveGroup> ZkSumcheckKey<G> {
    pub fn setup<H: HashToCurve<G>>(
        domain: &[u8],
        degree: usize,
        scalar: ScalarPedersen<G>,
    ) -> Result<Self, HashToCurveError> {
        assert!(degree >= 1, "round polynomial degree must be at least 1");

        let mut coeffs = VectorPedersen::setup::<H>(domain, degree + 1)?;
        coeffs.blind = scalar.blind;
        Ok(Self { coeffs, scalar })
    }

    /// この key が扱うラウンド多項式の次数上限。
    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ZkSumcheckProof<G: CurveGroup> {
    /// 各ラウンドの $s_i$ の係数への commitment $C_{s_i}$。
    pub poly_coms: Vec<G>,
    /// 各ラウンドの $s_i(r_i)$ への commitment $C_{e_i}$。
    pub eval_coms: Vec<G>,
    /// $\langle \mathrm{coeffs}(s_i), w_0 (2, 1, \ldots, 1) + w_1 (1, r_i, \ldots) \rangle
    /// = w_0 e_{i-1} + w_1 e_i$ の証明。
    pub proofs: Vec<LinearInnerProductProof<G>>,
}

/// `prove` の出力。proof に加えて、prover 側だけが持つ終点の情報を返す。
pub struct ZkSumcheckOutput<G: CurveGroup> {
    pub proof: ZkSumcheckProof<G>,
    /// 各ラウンドの challenge `r = (r_0, ..., r_{n-1})`。
    pub point: Vec<G::ScalarField>,
    /// 終点での各因子の値 `[p_0(r), ..., p_{k-1}(r)]`。
    pub final_evals: Vec<G::ScalarField>,
    /// 最終 claim $e_{n-1} = \mathrm{comb}(p_0(r), \ldots)$ と、その commitment の blind。
    pub final_claim: G::ScalarField,
    pub final_blind: G::ScalarField,
    /// 最終 claim の commitment。verifier が [`ZkSumcheckProof::verify`] で得るものと一致する。
    pub final_com: G,
}

/// 「$s(0) + s(1)$」と「$s(r)$」の二つの線形関数を $w_0, w_1$ で束ねた公開ベクトル。
fn combined_query<F: Field>(degree: usize, r: F, w0: F, w1: F) -> Vec<F> {
    let mut r_pow = F::one();
    (0..=degree)
        .map(|i| {
            let sum_weight = if i == 0 { F::from(2u64) } else { F::one() };
            let value = w0 * sum_weight + w1 * r_pow;
            r_pow *= r;
            value
        })
        .collect()
}

impl<G: CurveGroup> ZkSumcheckProof<G> {
    /// `claim` の commitment $\mathrm{Com}(\mathrm{claim}; \mathrm{claim\_blind})$ から始めて、
    /// $\sum_x \mathrm{comb}(p_0(x), \ldots, p_{k-1}(x)) = \mathrm{claim}$ を zero knowledge で示す。
    ///
    /// `comb` と各ラウンドの計算は平文版と同じで、次数は `key.degree()` に固定される。
    #[allow(clippy::too_many_arguments)]
    pub fn prove<
        P: MultilinearPoly<G::ScalarField>,
        Comb: Fn(&[G::ScalarField]) -> G::ScalarField,
    >(
        key: &ZkSumcheckKey<G>,
        mut polys: Vec<P>,
        mut claim: G::ScalarField,
        mut claim_blind: G::ScalarField,
        comb: Comb,
        transcript: &mut Transcript,
        rng: &mut (impl Rng + CryptoRng),
    ) -> ZkSumcheckOutput<G> {
        let degree = key.degree();
        let num_rounds = SumcheckProof::assert_shape(&polys, degree);
        let mut claim_com = key.scalar.commit(&claim, &claim_blind);

        let mut poly_coms = Vec::with_capacity(num_rounds);
        let mut eval_coms = Vec::with_capacity(num_rounds);
        let mut proofs = Vec::with_capacity(num_rounds);
        let mut r = Vec::with_capacity(num_rounds);

        for _ in 0..num_rounds {
            let coeffs = SumcheckProof::round_polynomial(&polys, claim, degree, &comb).to_coeffs();
            let poly_blind = G::ScalarField::rand(rng);
            let poly_com = key.coeffs.commit(coeffs.coeffs(), &poly_blind);
            transcript.append_serializable(b"zk-sumcheck-poly-commitment", &poly_com);

            let r_i = transcript.challenge_field::<G::ScalarField>(b"zk-sumcheck-challenge");
            for p in polys.iter_mut() {
                p.fold(r_i);
            }

            // 次ラウンドの claim e_i = s_i(r_i) は commitment だけを送る。
            let eval = coeffs.eval(r_i);
            let eval_blind = G::ScalarField::rand(rng);
            let eval_com = key.scalar.commit(&eval, &eval_blind);
            transcript.append_serializable(b"zk-sumcheck-claim-commitment", &claim_com);
            transcript.append_serializable(b"zk-sumcheck-eval-commitment", &eval_com);

            // s_i(0) + s_i(1) = e_{i-1} と s_i(r_i) = e_i を w で束ねる。
            let w0 = transcript.challenge_field::<G::ScalarField>(b"zk-sumcheck-weight");
            let w1 = transcript.challenge_field::<G::ScalarField>(b"zk-sumcheck-weight");
            let target = w0 * claim + w1 * eval;
            let target_blind = w0 * claim_blind + w1 * eval_blind;
            let target_com = claim_com * w0 + eval_com * w1;
            debug_assert_eq!(key.scalar.commit(&target, &target_blind), target_com);

            let query = combined_query(degree, r_i, w0, w1);
            let proof = LinearInnerProductProof::prove(
                &key.coeffs,
                &key.scalar,
                &poly_com,
                &target_com,
                coeffs.coeffs(),
                &poly_blind,
                &query,
                &target_blind,
                transcript,
                rng,
            );

            poly_coms.push(poly_com);
            eval_coms.push(eval_com);
            proofs.push(proof);
            r.push(r_i);
            claim = eval;
            claim_blind = eval_blind;
            claim_com = eval_com;
        }

        let final_evals = polys.iter().map(|p| p.final_constant()).collect();

        ZkSumcheckOutput {
            proof: Self {
                poly_coms,
                eval_coms,
                proofs,
            },
            point: r,
            final_evals,
            final_claim: claim,
            final_blind: claim_blind,
            final_com: claim_com,
        }
    }

    /// claim の commitment `claim_com` から始めて各ラウンドを検証し、
    /// 成功時は `(最終 claim の commitment, r)` を返す。
    ///
    /// ラウンド数は proof の長さから決まるので、呼び出し側は `r.len()` を確認すること。
    pub fn verify(
        &self,
        key: &ZkSumcheckKey<G>,
        claim_com: &G,
        transcript: &mut Transcript,
    ) -> Option<(G, Vec<G::ScalarField>)> {
        if self.poly_coms.len() != self.eval_coms.len() || self.poly_coms.len() != self.proofs.len()
        {
            return None;
        }

        let degree = key.degree();
        let mut claim_com = *claim_com;
        let mut r = Vec::with_capacity(self.poly_coms.len());

        for ((poly_com, eval_com), proof) in
            self.poly_coms.iter().zip(&self.eval_coms).zip(&self.proofs)
        {
            transcript.append_serializable(b"zk-sumcheck-poly-commitment", poly_com);
            let r_i = transcript.challenge_field::<G::ScalarField>(b"zk-sumcheck-challenge");

            transcript.append_serializable(b"zk-sumcheck-claim-commitment", &claim_com);
            transcript.append_serializable(b"zk-sumcheck-eval-commitment", eval_com);
            let w0 = transcript.challenge_field::<G::ScalarField>(b"zk-sumcheck-weight");
            let w1 = transcript.challenge_field::<G::ScalarField>(b"zk-sumcheck-weight");
            let target_com = claim_com * w0 + *eval_com * w1;

            let query = combined_query(degree, r_i, w0, w1);
            if !proof.verify(
                &key.coeffs,
                &key.scalar,
                &query,
                poly_com,
                &target_com,
                transcript,
            ) {
                return None;
            }

            r.push(r_i);
            claim_com = *eval_com;
        }

        Some((claim_com, r))
    }
}

#[cfg(test)]
mod tests {
    use super::{ZkSumcheckKey, ZkSumcheckOutput, ZkSumcheckProof};
    use crate::primitive::{DenseMultilinearPoly, MultilinearPoly, ScalarPedersen, Transcript};
    use ark_bls12_381::{g1::Config, Fr as F, G1Projective};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    fn key(degree: usize) -> ZkSumcheckKey<G1Projective> {
        let scalar = ScalarPedersen::setup::<G1Hasher>(b"zk-sumcheck-test").unwrap();
        ZkSumcheckKey::setup::<G1Hasher>(b"zk-sumcheck-test-coeffs", degree, scalar).unwrap()
    }

    fn quadratic_example() -> (Vec<DenseMultilinearPoly<F>>, F) {
        let a = DenseMultilinearPoly::new([2u64, 3, 5, 7, 11, 13, 17, 19].map(F::from).to_vec(), 3);
        let b =
            DenseMultilinearPoly::new([23u64, 29, 31, 37, 41, 43, 47, 53].map(F::from).to_vec(), 3);
        let claim = (0..8).map(|i| a.evals()[i] * b.evals()[i]).sum();
        (vec![a, b], claim)
    }

    #[test]
    fn zk_sumcheck_round_trips_and_hides_the_round_polynomials() {
        let key = key(2);
        let (polys, claim) = quadratic_example();
        let claim_blind = F::from(9);
        let claim_com = key.scalar.commit(&claim, &claim_blind);
        let comb = |v: &[F]| v[0] * v[1];

        let mut prover_transcript = Transcript::new(b"zk-sumcheck");
        let ZkSumcheckOutput {
            proof,
            point,
            final_evals,
            final_claim,
            final_blind,
            final_com,
        } = ZkSumcheckProof::prove(
            &key,
            polys.clone(),
            claim,
            claim_blind,
            comb,
            &mut prover_transcript,
            &mut StdRng::seed_from_u64(42),
        );

        assert_eq!(point.len(), 3);
        assert_eq!(
            final_evals,
            vec![polys[0].eval(&point), polys[1].eval(&point)]
        );
        assert_eq!(final_claim, comb(&final_evals));
        assert_eq!(final_com, key.scalar.commit(&final_claim, &final_blind));
        assert_eq!(&final_com, proof.eval_coms.last().unwrap());

        let mut verifier_transcript = Transcript::new(b"zk-sumcheck");
        let (verifier_final_com, verifier_point) = proof
            .verify(&key, &claim_com, &mut verifier_transcript)
            .unwrap();
        assert_eq!(verifier_point, point);
        assert_eq!(verifier_final_com, final_com);

        // 別の乱数で作った proof は commitment が全て異なるが、同じ statement を検証する。
        let mut other_transcript = Transcript::new(b"zk-sumcheck");
        let other = ZkSumcheckProof::prove(
            &key,
            polys,
            claim,
            claim_blind,
            comb,
            &mut other_transcript,
            &mut StdRng::seed_from_u64(7),
        );
        assert_ne!(other.proof.poly_coms, proof.poly_coms);
        let mut verifier_transcript = Transcript::new(b"zk-sumcheck");
        assert!(other
            .proof
            .verify(&key, &claim_com, &mut verifier_transcript)
            .is_some());
    }

    #[test]
    fn zero_rounds_return_the_initial_claim_commitment() {
        let key = key(2);
        let polys = vec![
            DenseMultilinearPoly::new(vec![F::from(3)], 0),
            DenseMultilinearPoly::new(vec![F::from(5)], 0),
        ];
        let (claim, claim_blind) = (F::from(15), F::from(4));
        let claim_com = key.scalar.commit(&claim, &claim_blind);

        let mut transcript = Transcript::new(b"zk-sumcheck");
        let output = ZkSumcheckProof::prove(
            &key,
            polys,
            claim,
            claim_blind,
            |v| v[0] * v[1],
            &mut transcript,
            &mut StdRng::seed_from_u64(42),
        );
        assert!(output.proof.poly_coms.is_empty());
        assert_eq!(output.final_com, claim_com);
        assert_eq!(output.final_blind, claim_blind);

        let mut transcript = Transcript::new(b"zk-sumcheck");
        let (final_com, point) = output
            .proof
            .verify(&key, &claim_com, &mut transcript)
            .unwrap();
        assert!(point.is_empty());
        assert_eq!(final_com, claim_com);
    }

    #[test]
    fn wrong_claim_commitment_and_tampered_proofs_are_rejected() {
        let key = key(2);
        let (polys, claim) = quadratic_example();
        let claim_blind = F::from(9);
        let claim_com = key.scalar.commit(&claim, &claim_blind);

        let mut transcript = Transcript::new(b"zk-sumcheck");
        let output = ZkSumcheckProof::prove(
            &key,
            polys,
            claim,
            claim_blind,
            |v| v[0] * v[1],
            &mut transcript,
            &mut StdRng::seed_from_u64(42),
        );

        // claim の値を変えた commitment。
        let wrong_claim_com = key.scalar.commit(&(claim + F::from(1)), &claim_blind);
        let mut transcript = Transcript::new(b"zk-sumcheck");
        assert!(output
            .proof
            .verify(&key, &wrong_claim_com, &mut transcript)
            .is_none());

        // 途中の eval commitment をずらす。
        let mut tampered = output.proof.clone();
        tampered.eval_coms[1] += key.scalar.generator;
        let mut transcript = Transcript::new(b"zk-sumcheck");
        assert!(tampered.verify(&key, &claim_com, &mut transcript).is_none());

        // 長さの整合しない proof。
        let mut tampered = output.proof.clone();
        tampered.proofs.pop();
        let mut transcript = Transcript::new(b"zk-sumcheck");
        assert!(tampered.verify(&key, &claim_com, &mut transcript).is_none());
    }

    #[test]
    fn key_degree_bounds_the_round_polynomials() {
        // 次数 3 の comb を次数 2 の key で証明すると、補間が壊れて自分の proof も通らない。
        let key = key(2);
        let polys = vec![DenseMultilinearPoly::new(
            [1u64, 2, 3, 4].map(F::from).to_vec(),
            2,
        )];
        let comb = |v: &[F]| v[0] * v[0] * v[0];
        let claim: F = (1..=4u64).map(|i| F::from(i * i * i)).sum();
        let claim_com = key.scalar.commit(&claim, &F::from(0));

        let mut transcript = Transcript::new(b"zk-sumcheck");
        let output = ZkSumcheckProof::prove(
            &key,
            polys,
            claim,
            F::from(0),
            comb,
            &mut transcript,
            &mut StdRng::seed_from_u64(42),
        );
        // 各ラウンドの内積関係自体は整合しているので verify は通るが、
        // 最終 claim が comb(final_evals) と一致しないので上位の最終検査で落ちる。
        let mut transcript = Transcript::new(b"zk-sumcheck");
        assert!(output
            .proof
            .verify(&key, &claim_com, &mut transcript)
            .is_some());
        assert_ne!(output.final_claim, comb(&output.final_evals));
    }
}
