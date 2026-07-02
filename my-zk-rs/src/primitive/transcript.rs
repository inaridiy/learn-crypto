use ark_ff::{Field, PrimeField};
use ark_serialize::CanonicalSerialize;
use sha2::{Digest, Sha256};

use super::matrix::Matrix;

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

    /// Arkworks の canonical serialization で field element を追加する。
    #[inline]
    pub fn append_field<F>(&mut self, label: &[u8], value: &F)
    where
        F: CanonicalSerialize,
    {
        self.append_serializable(label, value);
    }

    /// Arkworks の canonical serialization で値を transcript に追加する。
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

    /// Matrix の entry を row-major order で transcript に追加する。
    #[inline]
    pub fn append_matrix<F, const N: usize, const M: usize>(
        &mut self,
        label: &[u8],
        matrix: &Matrix<F, N, M>,
    ) where
        F: Field + CanonicalSerialize,
        [(); 1 << N]:,
        [(); 1 << M]:,
    {
        for row in matrix.rows() {
            for value in row {
                self.append_field(label, value);
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

    /// 現在の transcript から 32 byte の challenge を取り出し、その値を transcript に追加する。
    #[inline]
    pub fn challenge_bytes(&mut self, label: &[u8]) -> [u8; 32] {
        let digest = self.challenge_block(b"challenge-bytes", label, 0);
        self.append_bytes(label, &digest);
        digest
    }

    fn challenge_field_bytes<F: PrimeField>(&self, label: &[u8]) -> Vec<u8> {
        let byte_len = (F::MODULUS_BIT_SIZE as usize).div_ceil(8) + 16;
        let mut bytes = Vec::with_capacity(byte_len);
        let mut counter = 0u64;

        while bytes.len() < byte_len {
            let block = self.challenge_block(b"challenge-field", label, counter);
            bytes.extend_from_slice(&block);
            counter += 1;
        }

        bytes.truncate(byte_len);
        bytes
    }

    /// 現在の transcript から prime field element の challenge を取り出し、transcript に追加する。
    #[inline]
    pub fn challenge_field<F>(&mut self, label: &[u8]) -> F
    where
        F: PrimeField + CanonicalSerialize,
    {
        let bytes = self.challenge_field_bytes::<F>(label);
        let challenge = F::from_be_bytes_mod_order(&bytes);
        self.append_field(label, &challenge);
        challenge
    }
}

#[cfg(test)]
mod tests {
    use super::Transcript;

    use ark_bls12_381::Fr;

    #[test]
    fn same_transcript_produces_same_field_challenge() {
        let mut lhs = Transcript::new(b"example/v1");
        let mut rhs = Transcript::new(b"example/v1");

        lhs.append_bytes(b"message", b"hello");
        rhs.append_bytes(b"message", b"hello");

        assert_eq!(
            lhs.challenge_field::<Fr>(b"r"),
            rhs.challenge_field::<Fr>(b"r")
        );
    }

    #[test]
    fn different_messages_change_challenge() {
        let mut lhs = Transcript::new(b"example/v1");
        let mut rhs = Transcript::new(b"example/v1");

        lhs.append_bytes(b"message", b"hello");
        rhs.append_bytes(b"message", b"world");

        assert_ne!(
            lhs.challenge_field::<Fr>(b"r"),
            rhs.challenge_field::<Fr>(b"r")
        );
    }

    #[test]
    fn challenge_advances_transcript_state() {
        let mut transcript = Transcript::new(b"example/v1");

        let first = transcript.challenge_field::<Fr>(b"r");
        let second = transcript.challenge_field::<Fr>(b"r");

        assert_ne!(first, second);
    }
}
