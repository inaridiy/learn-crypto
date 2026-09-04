use ark_ff::{Field, PrimeField};
use ark_serialize::CanonicalSerialize;
use sha2::{Digest, Sha256};

use crate::primitive::Matrix;

/// Fiat-Shamir 変換で使う最小限の transcript。
///
/// Prover と verifier が同じ順序でメッセージを追加すると、同じ challenge を得る。
#[derive(Clone, Debug)]
pub struct Transcript {
    state: Sha256,
}

impl Transcript {
    /// Protocol ごとの domain separator で transcript を初期化する。
    #[inline]
    pub fn new(domain_separator: &[u8]) -> Self {
        let mut transcript = Self {
            state: Sha256::new(),
        };
        transcript.append_bytes(b"domain", domain_separator);
        transcript
    }

    /// ラベル付き byte 列を transcript に追加する。
    #[inline]
    pub fn append_bytes(&mut self, label: &[u8], bytes: &[u8]) {
        self.state.update((label.len() as u64).to_le_bytes());
        self.state.update(label);
        self.state.update((bytes.len() as u64).to_le_bytes());
        self.state.update(bytes);
    }

    /// `usize` を little-endian で transcript に追加する。
    #[inline]
    pub fn append_usize(&mut self, label: &[u8], value: usize) {
        self.append_bytes(label, &value.to_le_bytes());
    }

    #[inline]
    pub fn append_serializable<T>(&mut self, label: &[u8], value: &T)
    where
        T: CanonicalSerialize,
    {
        let mut bytes = Vec::with_capacity(value.uncompressed_size());
        value
            .serialize_uncompressed(&mut bytes)
            .expect("serializing into Vec cannot fail");
        self.append_bytes(label, &bytes);
    }

    /// 行列を形状ごと transcript に追加する。形状も bind しないと、
    /// 同じ成分列を持つ別形状の行列が同じ transcript になってしまう。
    #[inline]
    pub fn append_matrix<F, M>(&mut self, label: &[u8], matrix: &M)
    where
        F: Field,
        M: Matrix<F>,
    {
        self.append_usize(b"matrix-rows", matrix.rows());
        self.append_usize(b"matrix-cols", matrix.cols());
        for row in 0..matrix.rows() {
            for value in matrix.row(row) {
                self.append_serializable(label, &value);
            }
        }
    }

    fn challenge_block(&self, kind: &[u8], label: &[u8], counter: u64) -> [u8; 32] {
        let mut hasher = self.state.clone();
        hasher.update((kind.len() as u64).to_le_bytes());
        hasher.update(kind);
        hasher.update((label.len() as u64).to_le_bytes());
        hasher.update(label);
        hasher.update(counter.to_le_bytes());
        hasher.finalize().into()
    }

    /// 現在の transcript から challenge を導出し、その challenge 自体も transcript に追加する。
    #[inline]
    pub fn challenge_field<F>(&mut self, label: &[u8]) -> F
    where
        F: PrimeField + CanonicalSerialize,
    {
        // 法より 128 bit 長い入力を reduction し、単一の SHA-256 block を
        // そのまま field に写す場合の統計的偏りを negligible にする。
        let byte_len = (F::MODULUS_BIT_SIZE as usize).div_ceil(8) + 16;
        let mut bytes = Vec::with_capacity(byte_len);
        let mut counter = 0u64;

        while bytes.len() < byte_len {
            bytes.extend_from_slice(&self.challenge_block(b"challenge-field", label, counter));
            counter += 1;
        }
        bytes.truncate(byte_len);

        let challenge = F::from_be_bytes_mod_order(&bytes);
        self.append_serializable(label, &challenge);
        challenge
    }
}

#[cfg(test)]
mod tests {
    use super::Transcript;
    use ark_bls12_381::Fr;

    #[test]
    fn field_challenges_are_deterministic_and_advance_the_transcript() {
        let mut lhs = Transcript::new(b"transcript-test");
        let mut rhs = Transcript::new(b"transcript-test");
        lhs.append_bytes(b"message", b"hello");
        rhs.append_bytes(b"message", b"hello");

        let first = lhs.challenge_field::<Fr>(b"challenge");
        assert_eq!(first, rhs.challenge_field::<Fr>(b"challenge"));
        assert_ne!(first, lhs.challenge_field::<Fr>(b"challenge"));
    }

    #[test]
    fn field_challenges_bind_previous_messages() {
        let mut lhs = Transcript::new(b"transcript-test");
        let mut rhs = Transcript::new(b"transcript-test");
        lhs.append_bytes(b"message", b"hello");
        rhs.append_bytes(b"message", b"world");

        assert_ne!(
            lhs.challenge_field::<Fr>(b"challenge"),
            rhs.challenge_field::<Fr>(b"challenge")
        );
    }
}
