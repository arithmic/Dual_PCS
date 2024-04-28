#![allow(non_snake_case)]
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use crypto_bigint::Encoding;
use multilinear_kzg::common::{MleCommit, MleEvalProof};
use polynomial::LayerPolynomials;
use std::io::Write;

#[derive(Clone)]
pub struct GkrTranscript {
    pub final_evaluations: Vec<Vec<Scalar>>,
    pub claimed_values: Vec<Vec<MleLayerEvaluation>>,
    pub polynomials: Vec<LayerPolynomials>,
    pub final_layer_point: Vec<Scalar>,
}

impl GkrTranscript {
    pub fn new(
        final_evaluations: Vec<Vec<Scalar>>,
        claimed_values: Vec<Vec<MleLayerEvaluation>>,
        polynomials: Vec<LayerPolynomials>,
        final_layer_point: Vec<Scalar>,
    ) -> GkrTranscript {
        GkrTranscript {
            final_evaluations,
            claimed_values,
            polynomials,
            final_layer_point,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();

        //Write final evaluations
        for evaluations in &self.final_evaluations {
            for eval in evaluations {
                result
                    .write_all(&eval.0.to_be_bytes())
                    .expect("failed to write");
            }
        }

        //Write claimed_values
        for claimed_value in &self.claimed_values {
            for values in claimed_value {
                result
                    .write_all(&values.to_bytes())
                    .expect("failed to write");
            }
        }

        //Write polynomials
        for poly in &self.polynomials {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }
        //Note:- final layer points are not part of proof verifier can generate
        //final layer points using channel.
        result
    }
}
#[derive(Clone)]
pub struct CircuitLayer(Vec<Scalar>);

impl CircuitLayer {
    pub fn new(layer_evals: Vec<Scalar>) -> CircuitLayer {
        CircuitLayer(layer_evals)
    }

    pub fn get_position(&self, index: usize) -> Scalar {
        self.0[index]
    }
}

#[derive(Clone)]
pub struct CircuitBinaryTree {
    circuit_layers: Vec<CircuitLayer>,
    depth: usize,
}

impl CircuitBinaryTree {
    pub fn new(circuit_layers: Vec<CircuitLayer>) -> CircuitBinaryTree {
        let depth = circuit_layers.len() - 1;
        CircuitBinaryTree {
            circuit_layers,
            depth,
        }
    }
    pub fn get_layer_at_depth(&self, depth: usize) -> &CircuitLayer {
        &self.circuit_layers[depth]
    }
    pub fn get_depth(&self) -> usize {
        self.depth
    }
}

//Contains the evaluations of the MLE for the left and right child gate functions in the left and right position respectively.
#[derive(Clone)]
pub struct MleLayerEvaluation {
    pub left: Scalar,
    pub right: Scalar,
}

impl MleLayerEvaluation {
    pub fn new(left: Scalar, right: Scalar) -> MleLayerEvaluation {
        MleLayerEvaluation { left, right }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write left
        result
            .write_all(&self.left.0.to_be_bytes())
            .expect("failed to write");

        //Write right
        result
            .write_all(&self.right.0.to_be_bytes())
            .expect("failed to write");
        result
    }
}
#[derive(Clone)]

pub struct Commitments {
    pub A_commits: Vec<MleCommit<BlsCurve>>,
    pub B_commits: Vec<MleCommit<BlsCurve>>,
    pub C_commits: Vec<MleCommit<BlsCurve>>,
}

impl Commitments {
    pub fn new(
        A_commits: Vec<MleCommit<BlsCurve>>,
        B_commits: Vec<MleCommit<BlsCurve>>,
        C_commits: Vec<MleCommit<BlsCurve>>,
    ) -> Self {
        Self {
            A_commits,
            B_commits,
            C_commits,
        }
    }
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        result.push(mle_commit_to_bytes(&self.A_commits));
        result.push(mle_commit_to_bytes(&self.B_commits));
        result.push(mle_commit_to_bytes(&self.C_commits));
        result.len() as f64 / 1024f64
    }
}

pub fn mle_eval_proof_to_bytes(proofs: &Vec<MleEvalProof<BlsCurve>>) -> Vec<u8> {
    let mut result = Vec::new();
    for proof in proofs {
        result
            .write_all(&proof.evaluation.0.to_be_bytes())
            .expect("failed to write");
        for witness in &proof.witnesses {
            let affine_witness = witness.to_affine();
            result
                .write(&affine_witness.x.0.to_be_bytes())
                .expect("failed to write");
            result
                .write(&affine_witness.y.0.to_be_bytes())
                .expect("failed to write");
        }
    }
    result
}
pub fn mle_commit_to_bytes(commitments: &Vec<MleCommit<BlsCurve>>) -> Vec<u8> {
    let mut result = Vec::new();
    for commit in commitments {
        let affine_commit = commit.commitment.to_affine();
        result
            .write(&affine_commit.x.0.to_be_bytes())
            .expect("failed to write");
        result
            .write(&affine_commit.y.0.to_be_bytes())
            .expect("failed to write");
    }
    result
}
#[derive(Clone)]
pub struct Evaluations {
    pub A_evals: Vec<Scalar>,
    pub B_evals: Vec<Scalar>,
    pub C_evals: Vec<Scalar>,
}
impl Evaluations {
    pub fn new(A_evals: Vec<Scalar>, B_evals: Vec<Scalar>, C_evals: Vec<Scalar>) -> Self {
        Self {
            A_evals,
            B_evals,
            C_evals,
        }
    }
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        for value in &self.A_evals {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        for value in &self.B_evals {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        for value in &self.C_evals {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        result.len() as f64 / 1024f64
    }
}
