use bls381::scalar::Scalar;
use bls_curve::{
    bls::{AffinePoint, BlsCurve, ProjectivePoint},
    gt::Gt,
};
use crypto_bigint::Encoding;
use keccak::Keccak256;
use multilinear_kzg::common::MleCommit;
use non_algebraic_hash_traits::traits::NonAlgebraicHasher;
use random::RandomCoin;
use std::io::Write;
pub struct Channel {
    pub public_coin: RandomCoin<Keccak256<Scalar>>,
}

#[derive(Clone, Debug)]
pub struct MleEvalProofReseed {
    pub evaluation: Scalar,
    pub witnesses: Vec<AffinePoint>,
}

impl Channel {
    pub fn initialize_with_gt(public_inputs: &[Gt]) -> Self {
        let pub_inputs_bytes = public_inputs[0].0.to_bytes();
        let mut public_coin = RandomCoin::new(&pub_inputs_bytes);
        for item in public_inputs.iter().skip(1) {
            let final_bytes = item.0.to_bytes();
            let digest = Keccak256::<Scalar>::hash(&final_bytes);
            public_coin.reseed(digest);
        }
        Self { public_coin }
    }
    pub fn initialize_with_scalar(public_input: &[Scalar]) -> Self {
        let pub_inputs_bytes = Scalar::to_bytes(&public_input[0]);
        let mut public_coin = RandomCoin::new(&pub_inputs_bytes);
        for item in public_input.iter().skip(1) {
            let bytes = item.0.to_be_bytes();
            let digest = Keccak256::<Scalar>::hash(&bytes);
            public_coin.reseed(digest);
        }
        Self { public_coin }
    }
    pub fn initialize_with_affine_point(public_input: &[AffinePoint]) -> Channel {
        let pub_inputs_bytes = [
            public_input[0].x.0.to_be_bytes(),
            public_input[0].y.0.to_be_bytes(),
        ]
        .concat();
        let mut public_coin = RandomCoin::new(&pub_inputs_bytes);
        for item in public_input.iter().skip(1) {
            let bytes = [item.x.0.to_be_bytes(), item.y.0.to_be_bytes()].concat();
            let digest = Keccak256::<Scalar>::hash(&bytes);
            public_coin.reseed(digest);
        }
        Self { public_coin }
    }

    pub fn reseed_with_affine_point(&mut self, commits: &Vec<AffinePoint>) {
        commits.iter().for_each(|commit| {
            let mut bytes = Vec::new();
            bytes
                .write(&commit.x.0.to_be_bytes())
                .expect("failed to write");
            bytes
                .write(&commit.y.0.to_be_bytes())
                .expect("failed to write");
            let digest = Keccak256::<Scalar>::hash(&bytes);
            self.public_coin.reseed(digest);
        });
    }
    pub fn reseed_with_scalars(&mut self, data: &[Scalar]) {
        for item in data {
            let bytes = Scalar::to_bytes(item);
            let digest = Keccak256::<Scalar>::hash(&bytes);
            self.public_coin.reseed(digest);
        }
    }
    pub fn reseed_with_gt(&mut self, data: &[Gt]) {
        for item in data.iter() {
            let bytes = item.0.to_bytes();
            let digest = Keccak256::<Scalar>::hash(&bytes);
            self.public_coin.reseed(digest);
        }
    }

    pub fn get_random_point(&mut self) -> Scalar {
        self.public_coin.draw().unwrap()
    }

    pub fn get_random_points(&mut self, no_of_points: usize) -> Vec<Scalar> {
        let mut points = Vec::new();
        for _ in 0..no_of_points {
            points.push(self.get_random_point());
        }
        points
    }

    pub fn reseed_commits(&mut self, commits: Vec<Vec<ProjectivePoint>>) {
        for idx in 0..commits.len() {
            let final_bytes = bytes(&commits[idx]);
            let digest = Keccak256::<Scalar>::hash(&final_bytes);
            self.public_coin.reseed(digest);
        }
    }

    pub fn reseed_mle_commit(&mut self, mle_commit_vec: Vec<MleCommit<BlsCurve>>) {
        for idx in 0..mle_commit_vec.len() {
            let final_bytes = bytes(&vec![mle_commit_vec[idx].commitment.clone()]);
            let digest = Keccak256::<Scalar>::hash(&final_bytes);
            self.public_coin.reseed(digest);
        }
    }

    pub fn reseed_mle_eval_proof(&mut self, mle_eval_proof: &[MleEvalProofReseed]) {
        for idx in 0..mle_eval_proof.len() {
            let mut bytes = Vec::new();
            bytes
                .write(&mle_eval_proof[idx].evaluation.0.to_be_bytes())
                .expect("failed to write");
            for idx in 0..mle_eval_proof[idx].witnesses.len() {
                bytes
                    .write(&mle_eval_proof[idx].witnesses[idx].x.0.to_be_bytes())
                    .expect("failed to write");
                bytes
                    .write(&mle_eval_proof[idx].witnesses[idx].y.0.to_be_bytes())
                    .expect("failed to write");
            }
            let digest = Keccak256::<Scalar>::hash(&bytes);
            self.public_coin.reseed(digest);
        }
    }
}
pub fn bytes(inputs: &Vec<ProjectivePoint>) -> Vec<u8> {
    let mut final_bytes = Vec::new();
    for input in inputs {
        let affine_input = input.to_affine();
        let bytes = affine_input.x.0.to_be_bytes();
        for byte in bytes {
            final_bytes.push(byte);
        }
        let bytes = affine_input.y.0.to_be_bytes();
        for byte in bytes {
            final_bytes.push(byte);
        }
    }
    final_bytes
}
