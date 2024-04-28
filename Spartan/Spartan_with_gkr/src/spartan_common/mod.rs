#![allow(unused)]
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use crypto_bigint::{Encoding, Limb, Uint};
use grand_product_with_gkr::gkr_common::{CircuitBinaryTree, GkrTranscript};
use multilinear_kzg::common::{MleCommit, MleEvalProof};
use polynomial::{MultPolynomial, Polynomial};
use preprocessing::{mle_commit_to_bytes, mle_eval_proof_to_bytes, SparseCommit, SparseMetaData};
use std::io::Write;
use traits::traits::Field;

use crate::preprocessing;
extern crate bls381;
extern crate bls_curve;
extern crate crypto_bigint;
extern crate polynomial;
extern crate traits;
#[allow(unused)]
const SIX_INV: Scalar = Scalar(Uint {
    limbs: [
        Limb(15372286724512153601),
        Limb(1954008828163476649),
        Limb(9224930440102993583),
        Limb(6961264049553707793),
    ],
});
const TWO_INV: Scalar = Scalar(Uint {
    limbs: [
        Limb(9223372034707292161),
        Limb(12240451741123816959),
        Limb(1845609449319885826),
        Limb(4176758429732224676),
    ],
});

pub fn interpolate(eval: &mut Vec<Scalar>) {
    let t0 = eval[0] - eval[1].double();
    eval[1] = (-(eval[0] + eval[2]) - t0.double()) * TWO_INV;
    eval[2] = (t0 + eval[2]) * TWO_INV;
}
#[derive(Clone)]
#[allow(non_snake_case, non_camel_case_types)]
pub struct eR1CSmetadata {
    pub A: SparseMetaData,
    pub B: SparseMetaData,
    pub C: SparseMetaData,
}
#[allow(non_snake_case)]
impl eR1CSmetadata {
    pub fn new(A: SparseMetaData, B: SparseMetaData, C: SparseMetaData) -> eR1CSmetadata {
        eR1CSmetadata { A, B, C }
    }
}

#[allow(non_snake_case, non_camel_case_types)]
#[derive(Clone)]
pub struct eR1CSCommitments {
    pub A: SparseCommit,
    pub B: SparseCommit,
    pub C: SparseCommit,
    pub E: MleCommit<BlsCurve>,
    pub W: MleCommit<BlsCurve>,
}

#[allow(non_snake_case)]
impl eR1CSCommitments {
    pub fn new(
        A: SparseCommit,
        B: SparseCommit,
        C: SparseCommit,
        E: MleCommit<BlsCurve>,
        W: MleCommit<BlsCurve>,
    ) -> eR1CSCommitments {
        eR1CSCommitments { A, B, C, E, W }
    }
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        result
            .write_all(&self.A.to_bytes())
            .expect("failed to write");
        result
            .write_all(&self.B.to_bytes())
            .expect("failed to write");
        result
            .write_all(&self.C.to_bytes())
            .expect("failed to write");
        result
            .write_all(&mle_commit_to_bytes(&[self.E].to_vec()))
            .expect("failed to write");
        result
            .write_all(&mle_commit_to_bytes(&[self.W].to_vec()))
            .expect("failed to write");
        result.len() as f64 / 1024f64
    }
}
#[derive(Clone)]
#[allow(non_snake_case, non_camel_case_types)]
pub struct InitialSumCheckTranscript {
    pub polynomials: Vec<Polynomial>,
    pub random_points: Vec<Scalar>,
    pub Az_claimed_val: Scalar,
    pub Bz_claimed_val: Scalar,
    pub Cz_claimed_val: Scalar,
    pub E_claimed_val: Scalar,
}
#[allow(non_snake_case)]
impl InitialSumCheckTranscript {
    pub fn new(
        polynomials: Vec<Polynomial>,
        random_points: Vec<Scalar>,
        Az_claimed_val: Scalar,
        Bz_claimed_val: Scalar,
        Cz_claimed_val: Scalar,
        E_claimed_val: Scalar,
    ) -> InitialSumCheckTranscript {
        InitialSumCheckTranscript {
            polynomials,
            random_points,
            Az_claimed_val,
            Bz_claimed_val,
            Cz_claimed_val,
            E_claimed_val,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write poplynomials
        for poly in &self.polynomials {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write Az_claimed_val
        result
            .write_all(&self.Az_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write Bz_claimed_val
        result
            .write_all(&self.Bz_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write Cz_claimed_val
        result
            .write_all(&self.Cz_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write E_claimed_val
        result
            .write_all(&self.E_claimed_val.0.to_be_bytes())
            .expect("failed to write");
        result
    }
}

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct ParSumCheckTranscript {
    pub polynomials: Vec<Polynomial>,
    pub random_points: Vec<Scalar>,
    pub A_claimed_val: Scalar,
    pub B_claimed_val: Scalar,
    pub C_claimed_val: Scalar,
    pub Z_claimed_val: Scalar,
}

#[allow(non_snake_case)]
impl ParSumCheckTranscript {
    pub fn new(
        polynomials: Vec<Polynomial>,
        random_points: Vec<Scalar>,
        A_claimed_val: Scalar,
        B_claimed_val: Scalar,
        C_claimed_val: Scalar,
        Z_claimed_val: Scalar,
    ) -> ParSumCheckTranscript {
        ParSumCheckTranscript {
            polynomials,
            random_points,
            A_claimed_val,
            B_claimed_val,
            C_claimed_val,
            Z_claimed_val,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write poplynomials
        for poly in &self.polynomials {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write A_claimed_val
        result
            .write_all(&self.A_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write B_claimed_val
        result
            .write_all(&self.B_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write C_claimed_val
        result
            .write_all(&self.C_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write Z_claimed_val
        result
            .write_all(&self.Z_claimed_val.0.to_be_bytes())
            .expect("failed to write");
        result
    }
}

#[allow(non_snake_case, non_camel_case_types)]
#[derive(Clone)]
pub struct eR1CStranscript {
    pub first_sum_check_transcript: InitialSumCheckTranscript,
    pub par_sum_check_transcript: ParSumCheckTranscript,
    pub BatchProof: BatchSparseEvalProof,
    pub E_eval_proof: MleEvalProof<BlsCurve>,
    pub W_eval_proof: MleEvalProof<BlsCurve>,
}

#[allow(non_snake_case)]
impl eR1CStranscript {
    pub fn new(
        first_sum_check_transcript: InitialSumCheckTranscript,
        par_sum_check_transcript: ParSumCheckTranscript,
        BatchProof: BatchSparseEvalProof,
        E_eval_proof: MleEvalProof<BlsCurve>,
        W_eval_proof: MleEvalProof<BlsCurve>,
    ) -> eR1CStranscript {
        eR1CStranscript {
            first_sum_check_transcript,
            par_sum_check_transcript,
            BatchProof,
            E_eval_proof,
            W_eval_proof,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write first_sum_check_transcript
        result
            .write_all(&self.first_sum_check_transcript.to_bytes())
            .expect("failed to write");

        //Write par_sum_check_transcript
        result
            .write_all(&self.par_sum_check_transcript.to_bytes())
            .expect("failed to write");

        //Write BatchProof
        result
            .write_all(&self.BatchProof.to_bytes())
            .expect("failed to write");

        //Write E_eval_proof
        result
            .write_all(&mle_eval_proof_to_bytes(
                &[self.E_eval_proof.clone()].to_vec(),
            ))
            .expect("failed to write");

        //Write W_eval_proof
        result
            .write_all(&mle_eval_proof_to_bytes(
                &[self.W_eval_proof.clone()].to_vec(),
            ))
            .expect("failed to write");
        result
    }
}

#[derive(Clone)]
pub struct BatchSpartanSumCheckTranscript {
    pub polynomials: Vec<Polynomial>,
    pub random_points: Vec<Scalar>,
    pub e_rx: Vec<MultPolynomial>,
    pub e_ry: Vec<MultPolynomial>,
    pub val: Vec<MultPolynomial>,
}

impl BatchSpartanSumCheckTranscript {
    pub fn new(
        polynomials: Vec<Polynomial>,
        random_points: Vec<Scalar>,
        e_rx: Vec<MultPolynomial>,
        e_ry: Vec<MultPolynomial>,
        val: Vec<MultPolynomial>,
    ) -> BatchSpartanSumCheckTranscript {
        BatchSpartanSumCheckTranscript {
            polynomials,
            random_points,
            e_rx,
            e_ry,
            val,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write polynomials
        for poly in &self.polynomials {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write e_rx
        for poly in &self.e_rx {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write e_ry
        for poly in &self.e_ry {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write val
        for poly in &self.val {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }
        result
    }
}
#[derive(Clone)]
pub struct BatchSparseEvalProof {
    pub gkr_transcript1: GkrTranscript,
    pub gkr_transcript2: GkrTranscript,
    pub e_rx_evals_sum_check: Vec<Scalar>,
    pub e_ry_evals_sum_check: Vec<Scalar>,
    pub val_evals_sum_check: Vec<Scalar>,
    pub sum_check_eval_proof: MleEvalProof<BlsCurve>,
    pub e_rx_commits: Vec<MleCommit<BlsCurve>>,
    pub e_ry_commits: Vec<MleCommit<BlsCurve>>,
    pub final_ts_evals_row_mem_check: Vec<Scalar>,
    pub final_ts_evals_col_mem_check: Vec<Scalar>,
    pub row_evals_mem_check: Vec<Scalar>,
    pub col_evals_mem_check: Vec<Scalar>,
    pub read_ts_evals_row_mem_check: Vec<Scalar>,
    pub read_ts_evals_col_mem_check: Vec<Scalar>,
    pub e_rx_evals_mem_check: Vec<Scalar>,
    pub e_ry_evals_mem_check: Vec<Scalar>,
    pub gkr_batch_eval_proof1: MleEvalProof<BlsCurve>,
    pub gkr_batch_eval_proof2: MleEvalProof<BlsCurve>,
    pub sum_check_transcript: BatchSpartanSumCheckTranscript,
    pub circuit2_depth: usize,
}

impl BatchSparseEvalProof {
    pub fn new(
        sum_check_transcript: BatchSpartanSumCheckTranscript,
        gkr_transcript1: GkrTranscript,
        gkr_transcript2: GkrTranscript,
        e_rx_evals_sum_check: Vec<Scalar>,
        e_ry_evals_sum_check: Vec<Scalar>,
        val_evals_sum_check: Vec<Scalar>,
        sum_check_eval_proof: MleEvalProof<BlsCurve>,
        e_rx_commits: Vec<MleCommit<BlsCurve>>,
        e_ry_commits: Vec<MleCommit<BlsCurve>>,
        final_ts_evals_row_mem_check: Vec<Scalar>,
        final_ts_evals_col_mem_check: Vec<Scalar>,
        row_evals_mem_check: Vec<Scalar>,
        col_evals_mem_check: Vec<Scalar>,
        read_ts_evals_row_mem_check: Vec<Scalar>,
        read_ts_evals_col_mem_check: Vec<Scalar>,
        e_rx_evals_mem_check: Vec<Scalar>,
        e_ry_evals_mem_check: Vec<Scalar>,
        gkr_batch_eval_proof1: MleEvalProof<BlsCurve>,
        gkr_batch_eval_proof2: MleEvalProof<BlsCurve>,
        circuit2_depth: usize,
    ) -> BatchSparseEvalProof {
        Self {
            gkr_transcript1,
            gkr_transcript2,
            e_rx_evals_sum_check,
            e_ry_evals_sum_check,
            val_evals_sum_check,
            sum_check_eval_proof,
            e_rx_commits,
            e_ry_commits,
            final_ts_evals_row_mem_check,
            final_ts_evals_col_mem_check,
            row_evals_mem_check,
            col_evals_mem_check,
            read_ts_evals_row_mem_check,
            read_ts_evals_col_mem_check,
            e_rx_evals_mem_check,
            e_ry_evals_mem_check,
            gkr_batch_eval_proof1,
            gkr_batch_eval_proof2,
            sum_check_transcript,
            circuit2_depth,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write gkr transcript 1
        result
            .write_all(&self.gkr_transcript1.to_bytes())
            .expect("failed to write");

        //Write gkr transcript 2
        result
            .write_all(&self.gkr_transcript2.to_bytes())
            .expect("failed to write");

        //Write e_rx_evals_sum_check
        for value in &self.e_rx_evals_sum_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write e_rx_evals_sum_check
        for value in &self.e_ry_evals_sum_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write val_evals_sum_check
        for value in &self.val_evals_sum_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write sum_check_eval_proof
        result
            .write_all(&mle_eval_proof_to_bytes(
                &[self.sum_check_eval_proof.clone()].to_vec(),
            ))
            .expect("failed to write");

        //Write e_rx_commits
        result
            .write_all(&mle_commit_to_bytes(&self.e_rx_commits))
            .expect("failed to write");

        //Write e_ry_commits
        result
            .write_all(&mle_commit_to_bytes(&self.e_ry_commits))
            .expect("failed to write");

        //Write final_ts_evals_row_mem_check
        for value in &self.final_ts_evals_row_mem_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write final_ts_evals_col_mem_checks
        for value in &self.final_ts_evals_col_mem_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write row_evals_mem_check
        for value in &self.row_evals_mem_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write col_evals_mem_check
        for value in &self.col_evals_mem_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write read_ts_evals_row_mem_check
        for value in &self.read_ts_evals_row_mem_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write read_ts_evals_col_mem_check
        for value in &self.read_ts_evals_col_mem_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write e_rx_evals_mem_check
        for value in &self.e_rx_evals_mem_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write e_ry_evals_mem_check
        for value in &self.e_ry_evals_mem_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write gkr_batch_eval_proof1
        result
            .write_all(&mle_eval_proof_to_bytes(
                &[self.gkr_batch_eval_proof1.clone()].to_vec(),
            ))
            .expect("failed to write");

        //Write gkr_batch_eval_proof2
        result
            .write_all(&mle_eval_proof_to_bytes(
                &[self.gkr_batch_eval_proof2.clone()].to_vec(),
            ))
            .expect("failed to write");

        //Write sum_check_transcript
        result
            .write(&self.sum_check_transcript.to_bytes())
            .expect("failed to write");

        //Write circuit2_depth
        result
            .write_all(&self.circuit2_depth.to_be_bytes())
            .expect("failed tow write");
        result
    }
}
pub struct BatchGrandProductCircuits {
    pub w_init_circuit: Vec<CircuitBinaryTree>,
    pub s_circuit: Vec<CircuitBinaryTree>,
    pub w_update_circuit: Vec<CircuitBinaryTree>,
    pub r_circuit: Vec<CircuitBinaryTree>,
}

impl BatchGrandProductCircuits {
    pub fn new(
        w_init_circuit: Vec<CircuitBinaryTree>,
        s_circuit: Vec<CircuitBinaryTree>,
        w_update_circuit: Vec<CircuitBinaryTree>,
        r_circuit: Vec<CircuitBinaryTree>,
    ) -> BatchGrandProductCircuits {
        BatchGrandProductCircuits {
            w_init_circuit,
            s_circuit,
            w_update_circuit,
            r_circuit,
        }
    }
}
