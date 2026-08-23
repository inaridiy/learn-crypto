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

    #[inline]
    pub fn append_matrix<F, M>(&mut self, label: &[u8], matrix: &M)
    where
        F: CanonicalSerialize,
        M: Matrix<F>,
    {
        for row in 0..matrix.rows() {
            for value in matrix.row(row) {
                self.append_serializable(label, value);
            }
        }
    }
}
