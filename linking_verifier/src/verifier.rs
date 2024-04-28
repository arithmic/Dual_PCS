use dory::dory_verifier::dory_verifier;
use helper::{LinkingProof, Setup};

pub fn linking_verifier(linking_proof: LinkingProof, setup: Setup) {
    dory_verifier(linking_proof.instance1, setup.clone());
    dory_verifier(linking_proof.instance2, setup.clone());
}
