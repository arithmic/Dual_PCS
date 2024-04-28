/******************************************************************************
 * Copyright (c) 2022 FOLIUM LABS PRIVATE LIMITED. and its affiliates.        *
 ******************************************************************************/

#[cfg(test)]
mod tests {
    use crate::Keccak256;
    use bls381::scalar::Scalar;
    use hash_traits::traits::{ElementHasher, NonAlgebraicHasher};

    #[test]
    fn hash_padding() {
        let b1 = [1_u8, 2, 3];
        let b2 = [1_u8, 2, 3, 0];
        // adding a zero bytes at the end of a byte string should result in a different hash
        let r1 = Keccak256::<Scalar>::hash(&b1);
        let r2 = Keccak256::<Scalar>::hash(&b2);
        assert_ne!(r1, r2);
    }

    #[test]
    fn hash_elements_padding() {
        let e1: [Scalar; 2] = [Scalar::ONE, Scalar::ONE + Scalar::ONE];
        let e2 = [e1[0], e1[1], Scalar::ZERO];
        // adding a zero element at the end of a list of elements should result in a different hash
        let r1 = Keccak256::<Scalar>::hash_elements(&e1);
        let r2 = Keccak256::<Scalar>::hash_elements(&e2);
        assert_ne!(r1, r2);
    }
}
