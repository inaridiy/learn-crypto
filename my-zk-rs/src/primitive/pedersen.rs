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
