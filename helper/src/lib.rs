extern crate bls381;
extern crate bls_curve;
extern crate crypto_bigint;
extern crate curve_traits;
extern crate pairing;
extern crate rayon;
extern crate traits;

use bls381::{
    fp::{self, Fp},
    fp12::Fp12,
    fp2::Fp2,
    fp6::Fp6,
    scalar::Scalar,
};
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve, gt::Gt};
use crypto_bigint::{Limb, Uint};
use curve_traits::{AffinePoint, ProjectivePoint};
use std::io::Write;
pub mod helper;
pub const GT: Gt = Gt(Fp12 {
    c0: Fp6 {
        c0: Fp2 {
            c0: Fp(Uint {
                limbs: [
                    Limb(12198553741706481255),
                    Limb(11602998788307096840),
                    Limb(11880741070476681737),
                    Limb(14648777714741633045),
                    Limb(6208170935726264425),
                    Limb(700231262397062287),
                ],
            }),
            c1: Fp(Uint {
                limbs: [
                    Limb(15118301284587687014),
                    Limb(11847256713795813219),
                    Limb(11063887216559327677),
                    Limb(6663032656387263051),
                    Limb(14908799407724351085),
                    Limb(1158433630501275059),
                ],
            }),
        },
        c1: Fp2 {
            c0: Fp(Uint {
                limbs: [
                    Limb(17104075383120609249),
                    Limb(6403502104081293018),
                    Limb(12430162978881461862),
                    Limb(5363509998812569428),
                    Limb(8881104172992534896),
                    Limb(968017914908377644),
                ],
            }),
            c1: Fp(Uint {
                limbs: [
                    Limb(4703336871247272034),
                    Limb(6344427475015105450),
                    Limb(10791303394370929266),
                    Limb(2663631613512519645),
                    Limb(15562083668502839522),
                    Limb(1521007986910129339),
                ],
            }),
        },
        c2: Fp2 {
            c0: Fp(Uint {
                limbs: [
                    Limb(6268644363917632007),
                    Limb(6497416171251237069),
                    Limb(706774429876889517),
                    Limb(8151050653349760171),
                    Limb(1970874966534801385),
                    Limb(598693385913103615),
                ],
            }),
            c1: Fp(Uint {
                limbs: [
                    Limb(17699397697463256361),
                    Limb(8254155138368499300),
                    Limb(7526117805400273386),
                    Limb(5337588917649939020),
                    Limb(13934420908459106854),
                    Limb(1069332948893783833),
                ],
            }),
        },
    },
    c1: Fp6 {
        c0: Fp2 {
            c0: Fp(Uint {
                limbs: [
                    Limb(1788066812242076247),
                    Limb(5853993770372580872),
                    Limb(5622480066118394026),
                    Limb(15006590000516455541),
                    Limb(14101829159759870274),
                    Limb(907594725561505379),
                ],
            }),
            c1: Fp(Uint {
                limbs: [
                    Limb(9491433057446813527),
                    Limb(13198761130102191454),
                    Limb(13782543225768332073),
                    Limb(12098000010369278576),
                    Limb(4026588527201840938),
                    Limb(753738637186552040),
                ],
            }),
        },
        c1: Fp2 {
            c0: Fp(Uint {
                limbs: [
                    Limb(13066388097362179807),
                    Limb(14873679341531347264),
                    Limb(1074057887983349298),
                    Limb(11909737590233434784),
                    Limb(16631544038528799965),
                    Limb(1356966303050345188),
                ],
            }),
            c1: Fp(Uint {
                limbs: [
                    Limb(11967175573146563989),
                    Limb(969058183564605404),
                    Limb(7172478249463461573),
                    Limb(3722602066622870282),
                    Limb(7527977823389056994),
                    Limb(1685176929295848849),
                ],
            }),
        },
        c2: Fp2 {
            c0: Fp(Uint {
                limbs: [
                    Limb(16753073207637582756),
                    Limb(3813272685356338721),
                    Limb(8887603224038402),
                    Limb(8134721535505007323),
                    Limb(12443625139069431715),
                    Limb(149558483571461475),
                ],
            }),
            c1: Fp(Uint {
                limbs: [
                    Limb(7460508214269129587),
                    Limb(9030920352603187779),
                    Limb(13269539598400253413),
                    Limb(5466289635821009038),
                    Limb(13142064019075172745),
                    Limb(1191022261430813373),
                ],
            }),
        },
    },
});
#[derive(Clone)]
pub struct DegreeBoundSetup {
    pub tau_1: Vec<curve_traits::ProjectivePoint<BlsCurve>>,
    pub tau_2: Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    pub delta_1_left_commitment: Vec<Gt>,
    pub delta_1_right_commitment: Vec<Gt>,
    pub delta_2_right_commitment: Vec<Gt>,
    pub tau_ipp: Vec<Gt>,
    pub matrix_commitments: Vec<MatrixCommitment>,
}
impl DegreeBoundSetup {
    pub fn setup_for_specific_degree(&self, degree: usize) -> Setup {
        let log2_degree_bound = self.tau_1.len().trailing_zeros() as usize;
        let log2_degree = degree.trailing_zeros() as usize;
        Setup {
            tau_1: self.tau_1[0..degree].to_vec(),
            tau_2: self.tau_2[0..degree].to_vec(),
            delta_1_left_commitment: self.delta_1_left_commitment
                [log2_degree_bound - log2_degree..]
                .to_vec(),
            delta_1_right_commitment: self.delta_1_right_commitment
                [log2_degree_bound - log2_degree..]
                .to_vec(),
            delta_2_right_commitment: self.delta_2_right_commitment
                [log2_degree_bound - log2_degree..]
                .to_vec(),
            tau_ipp: self.tau_ipp[log2_degree_bound - log2_degree..].to_vec(),
            matrix_commitments: self.matrix_commitments[log2_degree - 1].clone(),
        }
    }
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        //Write tau_1
        for idx in 0..self.tau_1.len() {
            let affine_tau = self.tau_1[idx].to_affine();
            result
                .write_all(&fp::Fp::to_bytes(&affine_tau.x))
                .expect("failed to write bytes");
            result
                .write_all(&fp::Fp::to_bytes(&affine_tau.y))
                .expect("failed to write bytes");
        }
        //Write tau_2
        for idx in 0..self.tau_2.len() {
            let affine_tau = self.tau_2[idx].to_affine();
            result
                .write_all(&Fp2::<Fp>::to_bytes(&affine_tau.x))
                .expect("failed to write bytes");
            result
                .write_all(&Fp2::<Fp>::to_bytes(&affine_tau.y))
                .expect("failed to write bytes");
        }

        //Write delta_1_left_commitment
        for idx in 0..self.delta_1_left_commitment.len() {
            result
                .write_all(&self.delta_1_left_commitment[idx].0.to_bytes())
                .expect("failed to write bytes");
        }
        //Write delta_1_right_commitment
        for idx in 0..self.delta_1_right_commitment.len() {
            result
                .write_all(&self.delta_1_right_commitment[idx].0.to_bytes())
                .expect("failed to write bytes");
        }
        //Write delta_2_right_commitment
        for idx in 0..self.delta_2_right_commitment.len() {
            result
                .write_all(&self.delta_2_right_commitment[idx].0.to_bytes())
                .expect("failed to write bytes");
        }
        //Write tau_ipp
        for idx in 0..self.tau_ipp.len() {
            result
                .write_all(&self.tau_ipp[idx].0.to_bytes())
                .expect("failed to write bytes");
        }

        //Write matrix_commitments
        for idx in 0..self.matrix_commitments.len() {
            result
                .write_all(&self.matrix_commitments[idx].to_bytes())
                .expect("failed to write");
        }

        result.len() as f64 / 1024f64
    }
}

#[derive(Clone)]
pub struct Setup {
    pub tau_1: Vec<curve_traits::ProjectivePoint<BlsCurve>>,
    pub tau_2: Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    pub delta_1_left_commitment: Vec<Gt>,
    pub delta_1_right_commitment: Vec<Gt>,
    pub delta_2_right_commitment: Vec<Gt>,
    pub tau_ipp: Vec<Gt>,
    pub matrix_commitments: MatrixCommitment,
}

#[derive(Clone)]
pub struct MatrixCommitment {
    pub matrix_row_commits: Vec<ProjectivePoint<BlsCurve>>,
    pub matrix_commitment: Gt,
}

impl MatrixCommitment {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        for idx in 0..self.matrix_row_commits.len() {
            let affine_commit = self.matrix_row_commits[idx].to_affine();
            result
                .write_all(&fp::Fp::to_bytes(&affine_commit.x))
                .expect("failed to write bytes");
            result
                .write_all(&fp::Fp::to_bytes(&affine_commit.y))
                .expect("failed to write bytes");
        }
        result
            .write_all(&self.matrix_commitment.0.to_bytes())
            .expect("failed to write");
        result
    }
}
#[derive(Clone)]
pub struct EvaluationProof {
    pub dory_proof: DoryProof,
    pub evaluation: Scalar,
}
impl EvaluationProof {
    pub fn to_bytes(&self) -> f64 {
        self.dory_proof.to_bytes() + (Scalar::to_bytes(&self.evaluation).len() as f64 / 1024f64)
    }
}
#[derive(Clone)]
pub struct LinkingProof {
    pub instance1: DoryProof,
    pub instance2: DoryProof,
}

impl LinkingProof {
    pub fn to_bytes(&self) -> f64 {
        let instance1_proof_size = self.instance1.to_bytes();
        let instance2_proof_size = self.instance2.to_bytes();
        instance1_proof_size + instance2_proof_size
    }
}
#[derive(Clone)]
pub struct WitnessCommitment {
    pub commitment_to_witness: Gt,
    pub commitment_to_univariate: Gt,
    pub g2_power_poly: Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    pub univariate_polynomial: Vec<Scalar>,
    pub g2_power_witness: Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
}

#[derive(Clone)]
pub struct WitnessCommitmentOfDoryEvaluationProof {
    pub commitment_to_witness: Gt,
    pub commitment_to_univariate: Gt,
    pub g2_power_poly: Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    pub univariate_polynomial: Vec<Scalar>,
    pub g2_power_witness: Vec<ProjectivePoint<G2BlsCurve>>,
}

#[derive(Clone)]
pub struct DoryProof {
    pub intermediate_layer_commitments: Vec<Vec<Gt>>,
    pub intermediate_d: Vec<Vec<Gt>>,
    pub intermediate_cross_products: Vec<Vec<Gt>>,
    pub final_values: (AffinePoint<BlsCurve>, AffinePoint<G2BlsCurve>),
}
impl DoryProof {
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        //Write intermediate_layer_commitments
        for commitments in &self.intermediate_layer_commitments {
            for commit in commitments {
                result
                    .write_all(&commit.0.to_bytes())
                    .expect("failed to write bytes");
            }
        }
        //Write intermediate_d
        for values in &self.intermediate_d {
            for v in values {
                result
                    .write_all(&v.0.to_bytes())
                    .expect("failed to write bytes");
            }
        }
        //Write intermediate_cross_products
        for values in &self.intermediate_cross_products {
            for v in values {
                result
                    .write_all(&v.0.to_bytes())
                    .expect("failed to write bytes");
            }
        }
        result
            .write_all(&fp::Fp::to_bytes(&self.final_values.0.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.final_values.0.y))
            .expect("failed to write bytes");
        result
            .write_all(&Fp2::<Fp>::to_bytes(&self.final_values.1.x))
            .expect("failed to write bytes");
        result
            .write_all(&Fp2::<Fp>::to_bytes(&self.final_values.1.y))
            .expect("failed to write bytes");

        result.len() as f64 / 1024f64
    }
}
