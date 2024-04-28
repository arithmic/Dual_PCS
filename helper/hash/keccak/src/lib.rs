use std::{marker::PhantomData, slice};
extern crate hash_traits;
use bls381::scalar::Scalar;
use crypto_bigint::Encoding;
use hash_traits::traits::{ByteDigest, ElementHasher, NonAlgebraicHasher};
use sha3::Digest;
extern crate sha3;
use traits::traits::{Field, PrimeField};
use utils::ByteWriter;
extern crate traits;
extern crate utils;
pub mod test;
extern crate bls381;
extern crate crypto_bigint;

// keccak WITH 256-BIT OUTPUT
// ================================================================================================
/// Implementation of the [Hasher](super::Hasher) trait for Keccak256 hash function with 256-bit
/// output.
pub struct Keccak256<F: Field + PrimeField + 'static>(PhantomData<F>);
impl<F: Field + PrimeField + 'static> NonAlgebraicHasher for Keccak256<F> {
    type Digest = ByteDigest<32>;
    fn hash(bytes: &[u8]) -> Self::Digest {
        ByteDigest(sha3::Keccak256::digest(bytes).into())
    }
    fn merge(values: &[Self::Digest; 2]) -> Self::Digest {
        ByteDigest(sha3::Keccak256::digest(ByteDigest::digests_as_bytes(values)).into())
    }
    fn merge_with_int(seed: Self::Digest, value: u64) -> Self::Digest {
        let mut data = [0; 64];
        data[..32].copy_from_slice(&seed.0);
        data[32..].copy_from_slice(&Scalar::from(value).0.to_be_bytes());
        ByteDigest(sha3::Keccak256::digest(&data).into())
    }
}

impl<F: Field + PrimeField + 'static> ElementHasher<F> for Keccak256<F> {
    fn hash_elements(elements: &[F]) -> Self::Digest {
        if F::IS_CANONICAL {
            // when element's internal and canonical representations are the same, we can hash
            // element bytes directly
            let p = elements.as_ptr();
            let len = elements.len() * F::ELEMENT_BYTES;
            let bytes = unsafe { slice::from_raw_parts(p as *const u8, len) };
            ByteDigest(sha3::Keccak256::digest(bytes).into())
        } else {
            // when elements' internal and canonical representations differ, we need to serialize
            // them before hashing
            let mut hasher = ShaHasher::new();
            hasher.write(elements);
            ByteDigest(hasher.finalize())
        }
    }
}

// SHA HASHER
// ================================================================================================
/// Wrapper around SHA3 hasher to implement [ByteWriter] trait for it.
struct ShaHasher(sha3::Keccak256);
impl ShaHasher {
    pub fn new() -> Self {
        Self(sha3::Keccak256::new())
    }
    pub fn finalize(self) -> [u8; 32] {
        self.0.finalize().into()
    }
}
impl ByteWriter for ShaHasher {
    fn write_u8(&mut self, value: u8) {
        self.0.update(&[value]);
    }
    fn write_u8_slice(&mut self, values: &[u8]) {
        self.0.update(values);
    }
}
