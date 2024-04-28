use crossbeam::channel::unbounded;
use dory::dory_prover::dory_prover;
use helper::{LinkingProof, Setup, WitnessCommitment};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub fn linking_prover(witness_commitments: WitnessCommitment, setup: Setup) -> LinkingProof {
    let proof_inputs = vec![
        (
            setup.clone(),
            setup.matrix_commitments.matrix_row_commits.clone(),
            [
                witness_commitments.commitment_to_witness,
                setup.matrix_commitments.matrix_commitment,
                witness_commitments.commitment_to_univariate,
            ]
            .to_vec(),
            witness_commitments.g2_power_poly.clone(),
        ),
        (
            setup.clone(),
            setup.tau_1,
            [
                witness_commitments.commitment_to_witness,
                setup.tau_ipp[0],
                witness_commitments.commitment_to_witness,
            ]
            .to_vec(),
            witness_commitments.g2_power_witness,
        ),
    ];
    let (tx, rx) = unbounded();
    (0..2).into_par_iter().for_each_with(tx.clone(), |tx, idx| {
        let proof_input = proof_inputs[idx].clone();
        let proof = dory_prover(proof_input.0, proof_input.1, proof_input.2, proof_input.3);
        tx.send((idx, proof)).unwrap();
    });
    drop(tx);
    let mut instance1 = Vec::new();
    let mut instance2 = Vec::new();
    for (idx, proof) in rx {
        match idx {
            0 => instance1.push(proof),
            1 => instance2.push(proof),
            _ => {
                println!("Invalid index")
            }
        }
    }

    LinkingProof {
        instance1: instance1[0].clone(),
        instance2: instance2[0].clone(),
    }
}
