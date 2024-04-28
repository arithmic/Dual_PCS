#![allow(non_snake_case)]
use bls381::{fp, scalar::Scalar};
use bls_curve::bls::AffinePoint;
use channel::Channel;
use crypto_bigint::Encoding;
use std::io::Write;
use traits::traits::{Field, PrimeField};
pub const OFFSET: Scalar = Scalar::GENERATOR;
#[derive(Clone)]
pub struct GkrTranscript {
    oddframe: OODFrame,
    quotient_poly_commit_1: AffinePoint,
    quotient_poly_commit_2: AffinePoint,
    trace_commits: Vec<AffinePoint>,
    composition_poly_commit: AffinePoint,
    trace_length: usize,
    final_layers: Vec<Scalar>,
}

impl GkrTranscript {
    pub fn new(
        ood_frame: OODFrame,
        quotient_poly_commit_1: AffinePoint,
        quotient_poly_commit_2: AffinePoint,
        trace_commits: Vec<AffinePoint>,
        composition_poly_commit: AffinePoint,
        trace_length: usize,
        final_layers: Vec<Scalar>,
    ) -> Self {
        Self {
            oddframe: ood_frame,
            quotient_poly_commit_1,
            quotient_poly_commit_2,
            trace_commits,
            composition_poly_commit,
            trace_length,
            final_layers,
        }
    }
    pub fn get_odd_frame(&self) -> &OODFrame {
        &self.oddframe
    }
    pub fn get_quotient_poly_commits(&self) -> (AffinePoint, AffinePoint) {
        (self.quotient_poly_commit_1, self.quotient_poly_commit_2)
    }
    pub fn get_trace_commits(&self) -> &Vec<AffinePoint> {
        &self.trace_commits
    }
    pub fn composition_poly_commit(&self) -> AffinePoint {
        self.composition_poly_commit
    }
    pub fn get_trace_length(&self) -> usize {
        self.trace_length
    }
    pub fn get_final_layers(&self) -> &Vec<Scalar> {
        &self.final_layers
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //write ood frame
        result
            .write_all(&self.oddframe.to_bytes())
            .expect("failed to write");

        //Write quotient_poly_commit_1
        result
            .write_all(&fp::Fp::to_bytes(&self.quotient_poly_commit_1.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.quotient_poly_commit_1.y))
            .expect("failed to write bytes");

        //Write quotient_poly_commit_2
        result
            .write_all(&fp::Fp::to_bytes(&self.quotient_poly_commit_2.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.quotient_poly_commit_2.y))
            .expect("failed to write bytes");

        // Write trace commits
        for commit in &self.trace_commits {
            result
                .write_all(&fp::Fp::to_bytes(&commit.x))
                .expect("failed to write bytes");
            result
                .write_all(&fp::Fp::to_bytes(&commit.y))
                .expect("failed to write bytes");
        }

        //Write composition_poly_commit
        result
            .write_all(&fp::Fp::to_bytes(&self.composition_poly_commit.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.composition_poly_commit.y))
            .expect("failed to write bytes");

        //Write trace length
        result
            .write_all(&self.trace_length.to_le_bytes())
            .expect("failed to write");

        //Write final_layers_1
        for value in &self.final_layers {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("faield to write");
        }

        result
    }
}

pub struct LeafLayers {
    leaf_layer_1: Vec<Vec<Scalar>>,
    leaf_layer_2: Vec<Vec<Scalar>>,
}
impl LeafLayers {
    pub fn new(leaf_layer_1: Vec<Vec<Scalar>>, leaf_layer_2: Vec<Vec<Scalar>>) -> Self {
        Self {
            leaf_layer_1,
            leaf_layer_2,
        }
    }
    pub fn leaf_layer_1(&self) -> &Vec<Vec<Scalar>> {
        &self.leaf_layer_1
    }
    pub fn leaf_layer_2(&self) -> &Vec<Vec<Scalar>> {
        &self.leaf_layer_2
    }
}
#[derive(Clone)]
pub struct Commitments {
    A_commits: Vec<AffinePoint>,
    B_commits: Vec<AffinePoint>,
    C_commits: Vec<AffinePoint>,
}
impl Commitments {
    pub fn new(
        A_commits: Vec<AffinePoint>,
        B_commits: Vec<AffinePoint>,
        C_commits: Vec<AffinePoint>,
    ) -> Self {
        Self {
            A_commits,
            B_commits,
            C_commits,
        }
    }
    pub fn get_A_commits(&self) -> &Vec<AffinePoint> {
        &self.A_commits
    }
    pub fn get_B_commits(&self) -> &Vec<AffinePoint> {
        &self.B_commits
    }

    pub fn get_C_commits(&self) -> &Vec<AffinePoint> {
        &self.C_commits
    }
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        for commit in &self.A_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write");
        }
        for commit in &self.B_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write")
        }
        for commit in &self.C_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write")
        }
        result.len() as f64 / 1024f64
    }
}
#[derive(Clone)]
pub struct Evaluations {
    pub A_at_z: Vec<Scalar>,
    pub A_at_gz: Vec<Scalar>,
    pub B_at_z: Vec<Scalar>,
    pub B_at_gz: Vec<Scalar>,
    pub C_at_z: Vec<Scalar>,
    pub C_at_gz: Vec<Scalar>,
}
impl Evaluations {
    pub fn new(
        A_at_z: Vec<Scalar>,
        A_at_gz: Vec<Scalar>,
        B_at_z: Vec<Scalar>,
        B_at_gz: Vec<Scalar>,
        C_at_z: Vec<Scalar>,
        C_at_gz: Vec<Scalar>,
    ) -> Self {
        Self {
            A_at_z,
            A_at_gz,
            B_at_z,
            B_at_gz,
            C_at_z,
            C_at_gz,
        }
    }
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        for value in &self.A_at_z {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        for value in &self.A_at_gz {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        for value in &self.B_at_z {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        for value in &self.B_at_gz {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        for value in &self.C_at_z {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        for value in &self.C_at_gz {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        result.len() as f64 / 1024f64
    }
}
pub struct FinalLayers {
    final_layer_1: Vec<Scalar>,
    final_layer_2: Vec<Scalar>,
}
impl FinalLayers {
    pub fn new(final_layer_1: Vec<Scalar>, final_layer_2: Vec<Scalar>) -> Self {
        Self {
            final_layer_1,
            final_layer_2,
        }
    }
    pub fn final_layer_1(&self) -> &Vec<Scalar> {
        &self.final_layer_1
    }
    pub fn final_layer_2(&self) -> &Vec<Scalar> {
        &self.final_layer_2
    }
}
#[derive(Clone)]
pub struct OODFrame {
    current_frame: Vec<Scalar>,
    next_frame: Vec<Scalar>,
    composition_frame: Scalar,
}
impl OODFrame {
    pub fn new(
        current_frame: Vec<Scalar>,
        next_frame: Vec<Scalar>,
        composition_frame: Scalar,
    ) -> Self {
        Self {
            current_frame,
            next_frame,
            composition_frame,
        }
    }
    pub fn get_current_frame(&self) -> &Vec<Scalar> {
        &self.current_frame
    }
    pub fn get_next_frame(&self) -> &Vec<Scalar> {
        &self.next_frame
    }
    pub fn get_composition_frame(&self) -> Scalar {
        self.composition_frame
    }
    pub fn reseed_with_ood_frame(&self, channel: &mut Channel) {
        channel.reseed_with_scalars(&self.current_frame);
        channel.reseed_with_scalars(&self.next_frame);
        channel.reseed_with_scalars([self.composition_frame].as_ref());
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write current_frame
        for value in &self.current_frame {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write next frame
        for value in &self.next_frame {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write composition frame
        result
            .write_all(&self.composition_frame.0.to_be_bytes())
            .expect("failed to write");
        result
    }
}

pub fn evaluate_constraints(
    current: &Vec<Scalar>,
    next: &Vec<Scalar>,
    n_circuits: usize,
    final_layers_1: &Vec<Scalar>,
    evaluations: &mut Vec<Scalar>,
) {
    //Boundary Constraint
    (0..n_circuits).for_each(|idx| evaluations[idx] = current[2 * idx] - current[(2 * idx) + 1]);
    //Transition Constraints
    (0..n_circuits).for_each(|idx| {
        evaluations[n_circuits + idx] =
            (current[(2 * idx) + 1] * next[2 * idx]) - next[(2 * idx) + 1]
    });
    //Boundary Constraint
    (0..n_circuits).for_each(|idx| {
        evaluations[(2 * n_circuits) + idx] = current[(2 * idx) + 1] - final_layers_1[idx];
    });
}

pub fn merge_constraints(
    evaluations: &Vec<Scalar>,
    x: Scalar,
    constraint_comp_coeffs: Vec<Scalar>,
    trace_length: usize,
    g_trace: Scalar,
    n_circuits: usize,
) -> Scalar {
    let x_1 = x - Scalar::ONE;
    let x_g_trace_trace = x - g_trace.power_by([(trace_length - 1) as u64, 0, 0, 0]);
    let x_trace = (x.power_by([trace_length as u64, 0, 0, 0]) - Scalar::ONE)
        * x_g_trace_trace.invert().unwrap();
    let mut result1 = Scalar::ZERO;
    for idx in 0..n_circuits {
        let adjustment =
            constraint_comp_coeffs[idx] + (constraint_comp_coeffs[3 * n_circuits + idx] * x);
        result1 = result1 + (adjustment * evaluations[idx]);
    }
    result1 *= (x_1).invert().unwrap();

    let mut result2 = Scalar::ZERO;
    for idx in n_circuits..2 * n_circuits {
        let adjustment = constraint_comp_coeffs[idx] + constraint_comp_coeffs[3 * n_circuits + idx];
        result2 = result2 + (adjustment * evaluations[idx]);
    }
    result2 *= (x_trace).invert().unwrap();

    let mut result3 = Scalar::ZERO;
    for idx in 2 * n_circuits..3 * n_circuits {
        let adjustment = constraint_comp_coeffs[n_circuits + idx]
            + (constraint_comp_coeffs[(3 * n_circuits) + idx] * x);
        result3 = result3 + (adjustment * evaluations[idx]);
    }
    result3 *= (x_g_trace_trace).invert().unwrap();
    result1 + result2 + result3
}
