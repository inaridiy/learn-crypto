use ark_ec::{
    hashing::{HashToCurve, HashToCurveError},
    AffineRepr, CurveGroup,
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
    pub mul_bases: Vec<G::Affine>,
    pub blind: G,
}

impl<G: CurveGroup> VectorPedersen<G> {
    pub fn new(generators: Vec<G>, mul_bases: Vec<G::Affine>, blind: G) -> Self {
        Self {
            generators,
            mul_bases,
            blind,
        }
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
        let mul_bases = G::batch_convert_to_mul_base(&generators);
        let blind = hasher.hash(b"blind")?.into_group();

        Ok(Self::new(generators, mul_bases, blind))
    }

    pub fn len(&self) -> usize {
        self.generators.len()
    }

    pub fn is_empty(&self) -> bool {
        self.generators.is_empty()
    }

    /// $\mathrm{Com}(\vec{v}; r) = \sum_i v_i G_i + r H$
    #[inline]
    pub fn commit(&self, values: &[G::ScalarField], r: &G::ScalarField) -> G {
        assert!(
            values.len() <= self.generators.len(),
            "too many values for this vector commitment"
        );
        G::msm_unchecked(&self.mul_bases[..values.len()], values) + self.blind * r
    }
}
