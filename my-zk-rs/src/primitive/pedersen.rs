use ark_ec::{
    AffineRepr, CurveGroup,
    hashing::{HashToCurve, HashToCurveError},
};

#[derive(Clone, Debug)]
pub struct Pedersen<G: CurveGroup, const SIZE: usize> {
    pub generators: [G; SIZE],
    pub blind_generator: G,
}

impl<G: CurveGroup, const SIZE: usize> Pedersen<G, SIZE> {
    pub fn setup<H: HashToCurve<G>>(domain: &[u8]) -> Result<Self, HashToCurveError> {
        let hasher = H::new(domain)?;

        let blind_generator: G = hasher.hash(b"pedersen:blind")?.into_group();
        let generators = (0..SIZE)
            .map(|i| {
                Ok(hasher
                    .hash(format!("pedersen:generator#{i}").as_bytes())?
                    .into_group())
            })
            .collect::<Result<Vec<_>, HashToCurveError>>()?
            .try_into()
            .unwrap();

        Ok(Self {
            generators,
            blind_generator,
        })
    }

    pub fn commit(&self, values: [G::ScalarField; SIZE], r: G::ScalarField) -> G {
        let mut result = self.blind_generator * r;

        for i in 0..SIZE {
            result += self.generators[i] * values[i];
        }

        result
    }

    pub fn verify(&self, commitment: G, values: [G::ScalarField; SIZE], r: G::ScalarField) -> bool {
        let rhs = self.commit(values, r);

        commitment == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::{Fr as F, G1Projective, g1::Config};
    use ark_ec::hashing::{curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use sha2::Sha256;

    type G1Hasher =
        MapToCurveBasedHasher<G1Projective, DefaultFieldHasher<Sha256, 128>, WBMap<Config>>;

    #[test]
    fn commit_verifies_with_matching_values_and_blind() {
        let pedersen = Pedersen::<G1Projective, 2>::setup::<G1Hasher>(b"pedersen-test").unwrap();
        let values = [F::from(3), F::from(35)];
        let blind = F::from(9);

        let commitment = pedersen.commit(values, blind);

        assert!(pedersen.verify(commitment, values, blind));
        assert!(!pedersen.verify(commitment, [F::from(3), F::from(36)], blind));
        assert!(!pedersen.verify(commitment, values, F::from(10)));
    }
}
