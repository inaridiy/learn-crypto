use ark_ec::{
    AffineRepr, CurveGroup,
    hashing::{HashToCurve, HashToCurveError},
};

#[derive(Clone, Debug)]
pub struct ScalarPedersen<G: CurveGroup> {
    pub generator: G,
    pub blind: G,
}

impl<G: CurveGroup> ScalarPedersen<G> {
    pub fn new(generator: G, blind: G) -> Self {
        Self { generator, blind }
    }

    pub fn setup<H: HashToCurve<G>>(domain: &[u8]) -> Result<Self, HashToCurveError> {
        let hasher = H::new(domain)?;

        let generator = hasher.hash(b"generator")?.into_group();
        let blind = hasher.hash(b"blind")?.into_group();

        Ok(Self::new(generator, blind))
    }

    #[inline]
    pub fn commit(&self, value: &G::ScalarField, r: &G::ScalarField) -> G {
        self.generator * value + self.blind * r
    }
}

#[derive(Clone, Debug)]
pub struct VectorPedersen<G: CurveGroup> {
    pub generators: Vec<G>,
    pub blind: G,
}

impl<G: CurveGroup> VectorPedersen<G> {
    pub fn new(generators: Vec<G>, blind: G) -> Self {
        Self { generators, blind }
    }

    pub fn setup<H: HashToCurve<G>>(
        domain: &[u8],
        num_generators: usize,
    ) -> Result<Self, HashToCurveError> {
        let hasher = H::new(domain)?;

        let generators = (0..num_generators)
            .map(|i| {
                Ok(hasher
                    .hash(format!("generator#{i}").as_bytes())?
                    .into_group())
            })
            .collect::<Result<Vec<G>, HashToCurveError>>()?;
        let blind = hasher.hash(b"blind")?.into_group();

        Ok(Self::new(generators, blind))
    }

    /// $\mathrm{Com}(\vec{v}; r) = \sum_i v_i G_i + r H$
    #[inline]
    pub fn commit(&self, values: &[G::ScalarField], r: &G::ScalarField) -> G {
        self.generators
            .iter()
            .zip(values)
            .fold(G::zero(), |sum, (g, v)| sum + *g * v)
            + self.blind * r
    }
}
