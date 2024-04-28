#![allow(unused)]
use std::{ collections::HashMap, fmt::Display, io::Write};
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use crypto_bigint::Encoding;
use grand_product_with_gkr::gkr_common::CircuitBinaryTree;
use multilinear_kzg::{common::{MleCommit, MleEvalProof, SRS}, prover::commit};
use polynomial::MultPolynomial;
use rayon::iter::{  IntoParallelIterator, ParallelIterator};

#[derive(Clone, Copy, Debug)]
pub struct ColumnData {
    pub column: usize,
    pub value: Scalar,
}

impl ColumnData {
    pub fn new(
        column:usize,
        value:Scalar
    )->ColumnData{
        ColumnData{
            column,
            value,
        }
    }
}
#[allow(unused_assignments)]
pub fn get_tuples(
    sparse_representation:&SparseRep,
    n_cols: usize

)->(MultPolynomial, MultPolynomial, MultPolynomial){
    let sparsity = sparse_representation.sparsity().next_power_of_two();
    let rows = sparse_representation.fourcoeffs.len();
    let cols = n_cols;
    
    let tuples: Vec<(Scalar,Scalar,Scalar)>= (0..rows).into_par_iter().map(|row|{
        let entries = sparse_representation.fourcoeffs.get(&row).unwrap();
        entries.iter().map(|entry| 
        (Scalar::from(row as u64), Scalar::from(entry.column as u64), entry.value)
        ).collect::<Vec<(Scalar, Scalar, Scalar)>>()
    }
    ).flatten().collect();
    
    
    let row_col_max = if rows > cols {rows} else {cols};
    let mut row = vec![Scalar::from(row_col_max as u64); sparsity];
    let mut col = vec![Scalar::ZERO; sparsity];
    let mut val = vec![Scalar::ZERO; sparsity];

    (&mut row, &mut col, &mut val, tuples).into_par_iter().for_each(|(row, col, val, tuple)|{
    
    *row = tuple.0;
    *col = tuple.1;
    *val = tuple.2;
    
    }
    );

    (
        MultPolynomial::new(row),
        MultPolynomial::new(col),
        MultPolynomial::new(val),
    )
}

pub fn get_timestamps(
    row:&MultPolynomial,
    col:&MultPolynomial,
    val:&MultPolynomial,
    dim:usize,
    
)->TimeStamps{
    let length = val.len();
    let mut read_ts_row = vec![Scalar::ZERO; length];
    let mut read_ts_col = vec![Scalar::ZERO; length];

    let mut final_ts_row = vec![Scalar::ZERO; dim];
    let mut final_ts_col = vec![Scalar::ZERO; dim];

    for i in 0..length{
        read_ts_row[i] = final_ts_row[row.get_coeff(i).0.as_words()[0] as usize];
        final_ts_row[row.get_coeff(i).0.as_words()[0] as usize] += Scalar::ONE; 

        read_ts_col[i] = final_ts_col[col.get_coeff(i).0.as_words()[0] as usize];
        final_ts_col[col.get_coeff(i).0.as_words()[0] as usize] += Scalar::ONE; 
    }

    TimeStamps::new(
        MultPolynomial::new(read_ts_row),
        MultPolynomial::new(read_ts_col),
        MultPolynomial::new(final_ts_row),
        MultPolynomial::new(final_ts_col)
    )
    
}

#[derive(Clone)]
pub struct TimeStamps{
    pub read_ts_row: MultPolynomial,
    pub read_ts_col: MultPolynomial,
    pub final_ts_row: MultPolynomial,
    pub final_ts_col: MultPolynomial
}

impl TimeStamps {
    
    pub fn new(
        read_ts_row: MultPolynomial,
        read_ts_col: MultPolynomial,
        final_ts_row: MultPolynomial,
        final_ts_col: MultPolynomial
    )->TimeStamps{
        TimeStamps{
            read_ts_row,
            read_ts_col,
            final_ts_row,
            final_ts_col,
        }
    }
}


pub struct GrandProductCircuits {
    pub w_init_circuit: CircuitBinaryTree,
    pub s_circuit: CircuitBinaryTree,
    pub w_update_circuit: CircuitBinaryTree,
    pub r_circuit: CircuitBinaryTree,
}

impl GrandProductCircuits {
    pub fn new(
        w_init_circuit: CircuitBinaryTree,
        s_circuit: CircuitBinaryTree,
        w_update_circuit: CircuitBinaryTree,
        r_circuit: CircuitBinaryTree,
    ) -> GrandProductCircuits {
        GrandProductCircuits {
            w_init_circuit,
            s_circuit,
            w_update_circuit,
            r_circuit,
        }
    }
}
#[derive(Clone)]
pub struct SparseMetaData{
    pub row: MultPolynomial,
    pub col: MultPolynomial,
    pub val: MultPolynomial,
    pub timestamps: TimeStamps,
}

impl SparseMetaData {
    pub fn generate(
        sparse_representation:&SparseRep,
        n_cols: usize
    )->SparseMetaData{
        let (row, col, val) = get_tuples(sparse_representation, n_cols);
        let timestamps = get_timestamps(&row, &col, &val, sparse_representation.dim(n_cols));
        SparseMetaData{
            row,
            col,
            val,
            timestamps,
        }
    }
    pub fn commit(
        &self,
        srs: &SRS<BlsCurve>
    )->SparseCommit{
        SparseCommit::new(&self, srs)
    }
}

#[derive(Debug, Clone)]
pub struct SparseRep{
    pub fourcoeffs: HashMap<usize, Vec<ColumnData>>,
}

impl SparseRep{
    
    pub fn new(
        fourcoeffs: HashMap<usize, Vec<ColumnData>>,
    )->SparseRep{

        SparseRep { fourcoeffs }
    }

   pub fn dim(
        &self,
        n_cols: usize
    )->usize{
        let rows = self.fourcoeffs.len();
        let cols = n_cols;
        let row_col_max = if rows>cols {rows} else {cols};
        row_col_max.next_power_of_two()
    }
    pub fn sparsity(
        &self
    )->usize{
        let mut nonzero_entries = 0;
        for i in 0..self.fourcoeffs.len(){
            nonzero_entries +=self.fourcoeffs.get(&i).unwrap().len()
        }
        nonzero_entries
    }

    pub fn evaluate(
        self,
        basis_evals:&Vec<Scalar>,
        dim:usize
    )->Scalar{
        self.fourcoeffs.into_par_iter()
        //Iterate over keys in hashmap, which correspond to row indices of the matrix representation of the sparse polynomial.
        .flat_map(|(row_idx, row_entries)|

            //For each row, iterate over entries and multiply with corresponding lagrange basis element and flatten returned results
            //into a vector of scalars.
            row_entries.iter().map(|coldata|
            
            basis_evals[(row_idx<<dim) + coldata.column] * coldata.value
        
         ).collect::<Vec<Scalar>>() 
         //Sum over all returned values to get final evaluation value as an inner product.
        ).reduce(||Scalar::ZERO, |acc,e| acc + e)

    }

    pub fn get_metadata(
        &self,
        n_cols: usize
    )->SparseMetaData{
        SparseMetaData::generate(&self, n_cols)
    }
    pub fn bind_row_variable(
        &self,
        basis_evals:&Vec<Scalar>,
        n_cols: usize
    )->Vec<Scalar>{

        let mut result = vec![Scalar::ZERO; self.dim(n_cols)];

        self.fourcoeffs.iter().for_each(|(row_idx, row_entries)|{
            row_entries.iter().for_each(|coldata|
                result[coldata.column]+= basis_evals[*row_idx]*coldata.value
            )

        });

        result
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SparseCommit{
    pub row_commit:      MleCommit<BlsCurve>,
    pub col_commit:         MleCommit<BlsCurve>,
    pub val_commit:         MleCommit<BlsCurve>,
    pub read_ts_row_commit: MleCommit<BlsCurve>,
    pub final_ts_row_commit: MleCommit<BlsCurve>,
    pub read_ts_col_commit: MleCommit<BlsCurve>,
    pub final_ts_col_commit: MleCommit<BlsCurve>,
}

impl SparseCommit {
    pub fn new(
       sparse_metadata:&SparseMetaData,
       srs: &SRS<BlsCurve>
    )->SparseCommit{
        let row_commit = commit(sparse_metadata.row.as_coeffs(), srs);                        
        let col_commit = commit(sparse_metadata.col.as_coeffs(), srs);                        
        let val_commit = commit(sparse_metadata.val.as_coeffs(), srs);                        
        let read_ts_row_commit   = commit(sparse_metadata.timestamps.read_ts_row.as_coeffs(),srs);
        let final_ts_row_commit = commit(sparse_metadata.timestamps.final_ts_row.as_coeffs(),srs);
        let read_ts_col_commit   = commit(sparse_metadata.timestamps.read_ts_col.as_coeffs(),srs);
        let final_ts_col_commit = commit(sparse_metadata.timestamps.final_ts_col.as_coeffs(),srs);
        
            SparseCommit{
                row_commit,
                col_commit,
                val_commit,
                read_ts_row_commit,
                final_ts_row_commit,
                read_ts_col_commit,
                final_ts_col_commit,
            } 
    }
    pub fn to_bytes(&self) -> Vec<u8>{
        let mut result = Vec::new();
        result.write_all(&mle_commit_to_bytes(&[self.row_commit].to_vec())).expect("failed to write");
        result.write_all(&mle_commit_to_bytes(&[self.col_commit].to_vec())).expect("failed to write");
        result.write_all(&mle_commit_to_bytes(&[self.val_commit].to_vec())).expect("failed to write");
        result.write_all(&mle_commit_to_bytes(&[self.read_ts_row_commit].to_vec())).expect("failed to write");
        result.write_all(&mle_commit_to_bytes(&[self.read_ts_col_commit].to_vec())).expect("failed to write");
        result.write_all(&mle_commit_to_bytes(&[self.final_ts_row_commit].to_vec())).expect("failed to write");
        result.write_all(&mle_commit_to_bytes(&[self.final_ts_col_commit].to_vec())).expect("failed to write");
        result
    }
}

impl Display for SparseCommit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} {} {} {} {} {} {}",
            self.row_commit, self.col_commit, self.val_commit, self.read_ts_row_commit,  self.final_ts_row_commit, self.read_ts_col_commit, self.final_ts_col_commit, 
        )
        
    }
}


impl SparseCommit{
    pub fn from_string(str: &str) -> Option<SparseCommit> {
     
        let mut split = str.split(" ");
        let row_commit_x = split.next()?;
        let row_commit_y = split.next()?;
        let row_commit_z = split.next()?;
        
        let row_commit_string = format!("{row_commit_x} {row_commit_y} {row_commit_z}");
        

        let col_commit_x = split.next()?;
        let col_commit_y = split.next()?;
        let col_commit_z = split.next()?;
        
        let col_commit_string = format!("{col_commit_x} {col_commit_y} {col_commit_z}");
        
        let val_commit_x = split.next()?;
        let val_commit_y = split.next()?;
        let val_commit_z = split.next()?;
        
        let val_commit_string = format!("{val_commit_x} {val_commit_y} {val_commit_z}");
       

        let read_ts_row_commit_x = split.next()?;
        let read_ts_row_commit_y = split.next()?;
        let read_ts_row_commit_z = split.next()?;
        
        let read_ts_row_commit_string = format!("{read_ts_row_commit_x} {read_ts_row_commit_y} {read_ts_row_commit_z}");
        

        let final_ts_row_commit_x = split.next()?;
        let final_ts_row_commit_y = split.next()?;
        let final_ts_row_commit_z = split.next()?;
        
        let final_ts_row_commit_string = format!("{final_ts_row_commit_x} {final_ts_row_commit_y} {final_ts_row_commit_z}");
        

        let read_ts_col_commit_x = split.next()?;
        let read_ts_col_commit_y = split.next()?;
        let read_ts_col_commit_z = split.next()?;
        
        let read_ts_col_commit_string = format!("{read_ts_col_commit_x} {read_ts_col_commit_y} {read_ts_col_commit_z}");
       

        let final_ts_col_commit_x = split.next()?;
        let final_ts_col_commit_y = split.next()?;
        let final_ts_col_commit_z = split.next()?;
        
        let final_ts_col_commit_string = format!("{final_ts_col_commit_x} {final_ts_col_commit_y} {final_ts_col_commit_z}");
        

        
        Some(
            SparseCommit{
                row_commit:  MleCommit::<BlsCurve>::from_string(&row_commit_string).unwrap() ,
                col_commit:  MleCommit::<BlsCurve>::from_string(&col_commit_string).unwrap(),
                val_commit:  MleCommit::<BlsCurve>::from_string(&val_commit_string).unwrap(),
                read_ts_row_commit:  MleCommit::<BlsCurve>::from_string(&read_ts_row_commit_string).unwrap(),
                final_ts_row_commit: MleCommit::<BlsCurve>::from_string(&final_ts_row_commit_string).unwrap(),
                read_ts_col_commit:  MleCommit::<BlsCurve>::from_string(&read_ts_col_commit_string).unwrap(),
                final_ts_col_commit:  MleCommit::<BlsCurve>::from_string(&final_ts_col_commit_string).unwrap(),
            }
        )
    }


 }

pub fn compute_coeff(r: &Vec<Scalar>) -> Vec<Scalar> {
    //Initialize fc_eq with (1- r[0]) and r[0]
    let mut fc = [Scalar::ONE - r[0], r[0]].to_vec();
    //Iterate over the length of the r vector
    for k in 1..r.len() {
        let temp = fc;
        fc = vec![Scalar::ZERO; temp.len() * 2];
        for iter in 0..temp.len() {
            fc[2 * iter] = temp[iter] * (Scalar::ONE - r[k as usize]);
            fc[2 * iter + 1] = temp[iter] * r[k as usize];
        }
    }
    fc
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

pub fn mle_eval_proof_to_bytes(proofs: &Vec<MleEvalProof<BlsCurve>>)-> Vec<u8> {
    let mut result = Vec::new();
    for proof in proofs{
        result.write_all(&proof.evaluation.0.to_be_bytes()).expect("failed to write");
       for witness in &proof.witnesses{
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