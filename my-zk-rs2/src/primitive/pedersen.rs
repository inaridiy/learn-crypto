use ark_ec::{
    AffineRepr, CurveGroup,
    hashing::{HashToCurve, HashToCurveError},
};

use super::helpers::msm_with_bases;

/// Pedersen vector commitment の生成元一式。
///
/// $\mathrm{Com}(\vec{v}; r) = \sum_i v_i G_i + r H$
///
#[derive(Clone, Debug)]
pub struct Pedersen<G: CurveGroup> {
    pub generators: Vec<G>,
    pub mul_bases: Vec<G::Affine>,
    pub blind: G,
}

impl<G: CurveGroup> Pedersen<G> {
    pub fn setup<H: HashToCurve<G>>(
        domain: &[u8],
        label: &str,
        num_generators: usize,
    ) -> Result<Self, HashToCurveError> {
        let hasher = H::new(domain)?;

        let generators = (0..num_generators)
            .map(|i| {
                Ok(hasher
                    .hash(format!("pedersen:{label}:g#{i}").as_bytes())?
                    .into_group())
            })
            .collect::<Result<Vec<G>, HashToCurveError>>()?;
        let mul_bases = G::batch_convert_to_mul_base(&generators);
        let blind = hasher.hash(b"pedersen:h")?.into_group();

        Ok(Self {
            generators,
            mul_bases,
            blind,
        })
    }

    /// $\mathrm{Com}(\vec{v}; r) = \sum_i v_i G_i + r H$
    #[inline]
    pub fn commit(&self, values: &[G::ScalarField], blind: &G::ScalarField) -> G {
        assert!(
            values.len() <= self.generators.len(),
            "too many values for this committer"
        );
        msm_with_bases::<G>(values, &self.mul_bases[..values.len()]) + self.blind * blind
    }

    /// スカラー 1 個のコミットメント $\mathrm{Com}(v; r) = v G_0 + r H$。
    #[inline]
    pub fn commit_scalar(&self, value: &G::ScalarField, blind: &G::ScalarField) -> G {
        self.commit(std::slice::from_ref(value), blind)
    }
}

#[cfg(test)]
mod tests {
    use super::Pedersen;
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    fn committer(label: &str, n: usize) -> Pedersen<G1Projective> {
        Pedersen::setup::<G1Hasher>(b"pedersen-test", label, n).unwrap()
    }

    #[test]
    fn commitment_is_additively_homomorphic() {
        let committer = committer("vec", 3);

        let lhs = committer.commit(&[F::from(1), F::from(2), F::from(3)], &F::from(7));
        let rhs = committer.commit(&[F::from(10), F::from(20), F::from(30)], &F::from(11));
        let sum = committer.commit(&[F::from(11), F::from(22), F::from(33)], &F::from(18));

        assert_eq!(lhs + rhs, sum);
    }

    #[test]
    fn committers_with_the_same_domain_share_the_blind_generator() {
        let lhs = committer("vec", 2);
        let rhs = committer("scalar", 1);

        assert_eq!(lhs.blind, rhs.blind);
        assert_ne!(lhs.generators[0], rhs.generators[0]);
    }

    #[test]
    fn commit_scalar_uses_the_first_generator() {
        let committer = committer("scalar", 1);
        let (value, blind) = (F::from(5), F::from(7));

        assert_eq!(
            committer.commit_scalar(&value, &blind),
            committer.generators[0] * value + committer.blind * blind
        );
    }
}
