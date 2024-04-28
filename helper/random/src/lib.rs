pub mod errors;
use bls381::{fp::Fp, scalar::Scalar};
use errors::RandomCoinError;
use hash_traits::traits::{Digest, NonAlgebraicHasher};
use traits::traits::PrimeField;

use core::{convert::TryInto, marker::PhantomData};
use crypto_bigint::U256;
use utils::collections::Vec;
pub const SCALAR_MODULUS: U256 = U256::from_be_hex(Scalar::MODULUS);

// RANDOM COIN
// ================================================================================================
/// Pseudo-random element generator for finite fields.
///
/// A random coin can be used to draws elements uniformly at random from the specified base field
// (which is specified via the `B` type parameter) or from any extension of the base field.
///
/// Internally we use a cryptographic hash function (which is specified via the `H` type parameter),
/// to draw elements from the field. The coin works roughly as follows:
/// - The internal state of the coin consists of a `seed` and a `counter`. At instantiation
///   time, the `seed` is set to a hash of the provided bytes, and the `counter` is set to 0.
/// - To draw the next element, we increment the `counter` and compute hash(`seed` || `counter`).
///   If the resulting value is a valid field element, we return the result; otherwise we try
///   again until a valid element is found or the number of allowed tries is exceeded.
/// - We can also re-seed the coin with a new value. During the reseeding procedure, the
///   seed is set to hash(`old_seed` || `new_seed`), and the counter is reset to 0.
///
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RandomCoin<H>
where
    H: NonAlgebraicHasher,
{
    pub seed: H::Digest,
    pub counter: u64,
    _base_field: PhantomData<Fp>,
}
impl<H: NonAlgebraicHasher> RandomCoin<H> {
    // CONSTRUCTOR
    // --------------------------------------------------------------------------------------------
    /// Returns a new random coin instantiated with the provided `seed`.
    pub fn new(seed: &[u8]) -> Self {
        let seed = H::hash(seed);
        RandomCoin {
            seed,
            counter: 0,
            _base_field: PhantomData,
        }
    }
    // RESEEDING
    // --------------------------------------------------------------------------------------------
    /// Reseeds the coin with the specified data by setting the new seed to hash(`seed` || `data`).
    ///
    /// # Examples
    /// ```
    pub fn reseed(&mut self, data: H::Digest) {
        self.seed = H::merge(&[self.seed, data]);
        self.counter = 0;
    }
    /// Reseeds the coin with the specified value by setting the new seed to hash(`seed` ||
    /// `value`).
    ///
    /// # Examples
    /// ```
    /// ```
    pub fn reseed_with_int(&mut self, value: u64) {
        self.seed = H::merge_with_int(self.seed, value);
        self.counter = 0;
    }
    // PUBLIC ACCESSORS
    // --------------------------------------------------------------------------------------------
    /// Returns the number of leading zeros in the seed if it is interpreted as an integer in
    /// big-endian byte order.
    ///
    /// # Examples
    /// ```
    /// ```
    pub fn leading_zeros(&self) -> u32 {
        let bytes = self.seed.as_bytes();
        let seed_head = u64::from_le_bytes(bytes[..8].try_into().unwrap());
        seed_head.trailing_zeros()
    }
    /// Computes hash(`seed` || `value`) and returns the number of leading zeros in the resulting
    /// value if it is interpreted as an integer in big-endian byte order.
    pub fn check_leading_zeros(&self, value: u64) -> u32 {
        let new_seed = H::merge_with_int(self.seed, value);
        let bytes = new_seed.as_bytes();
        let seed_head = u64::from_le_bytes(bytes[..8].try_into().unwrap());
        seed_head.trailing_zeros()
    }
    // DRAW METHODS
    // --------------------------------------------------------------------------------------------
    /// Returns the next pseudo-random field element.
    ///
    /// # Errors
    /// Returns an error if a valid field element could not be generated after 1000 calls to the
    /// PRNG.
    pub fn draw(&mut self) -> Result<Scalar, RandomCoinError>
where {
        for _ in 0..1000 {
            // get the next pseudo-random value and take the first ELEMENT_BYTES from it
            let value = self.next();
            let bytes = &value.as_bytes()[..std::mem::size_of::<U256>()];
            let uint = U256::from_be_slice(bytes.into());
            let element = Scalar(uint >> 3);
            return Ok(element);
        }
        Err(RandomCoinError::FailedToDrawFieldElement(1000))
    }
    /// Returns the next pair of pseudo-random field elements.
    ///
    /// # Errors
    /// Returns an error if any of the field elements could not be generated after 100 calls to
    /// the PRNG;
    pub fn draw_pair(&mut self) -> Result<(Scalar, Scalar), RandomCoinError> {
        Ok((self.draw()?, self.draw()?))
    }
    /// Returns the next triplet of pseudo-random field elements.
    ///
    /// # Errors
    /// Returns an error if any of the field elements could not be generated after 100 calls to
    /// the PRNG;
    pub fn draw_triple(&mut self) -> Result<(Scalar, Scalar, Scalar), RandomCoinError> {
        Ok((self.draw()?, self.draw()?, self.draw()?))
    }
    /// Returns a vector of unique integers selected from the range [0, domain_size).
    ///
    /// # Errors
    /// Returns an error if the specified number of unique integers could not be generated
    /// after 1000 calls to the PRNG.
    ///
    /// # Panics
    /// Panics if:
    /// - `domain_size` is not a power of two.
    /// - `num_values` is greater than or equal to `domain_size`.
    ///
    pub fn draw_integers(
        &mut self,
        num_values: usize,
        domain_size: usize,
    ) -> Result<Vec<usize>, RandomCoinError> {
        assert!(
            domain_size.is_power_of_two(),
            "domain size must be a power of two"
        );
        assert!(
            num_values < domain_size,
            "number of values must be smaller than domain size"
        );
        // determine how many bits are needed to represent valid values in the domain
        let v_mask = (domain_size - 1) as u64;
        // draw values from PRNG until we get as many unique values as specified by num_queries
        let mut values = Vec::new();
        for _ in 0..1000 {
            // get the next pseudo-random value and read the first 8 bytes from it
            let bytes: [u8; 8] = self.next().as_bytes()[..8].try_into().unwrap();
            // convert to integer and limit the integer to the number of bits which can fit
            // into the specified domain
            let value = (u64::from_le_bytes(bytes) & v_mask) as usize;
            if values.contains(&value) {
                continue;
            }
            values.push(value);
            if values.len() == num_values {
                break;
            }
        }
        if values.len() < num_values {
            return Err(RandomCoinError::FailedToDrawIntegers(
                num_values,
                values.len(),
                1000,
            ));
        }
        Ok(values)
    }
    // HELPER METHODS
    // --------------------------------------------------------------------------------------------
    /// Updates the state by incrementing the counter and returns hash(seed || counter)
    fn next(&mut self) -> H::Digest {
        self.counter += 1;
        H::merge_with_int(self.seed, self.counter)
    }
}
