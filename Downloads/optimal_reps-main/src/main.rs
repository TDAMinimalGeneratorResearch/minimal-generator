use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind};
use exhact::clique::Simplex;
use num::rational::Ratio;
use std;
use coin_cbc::{raw::Status, Col, Model, Sense, Solution as CbcSolution};
use sprs::{CsMat, CsVec};
use std::collections::HashMap;



fn main() {
    // Create the problem.
    let mut m = Model::default();
    m.set_parameter("log", "0"); // turn off logging 
    let obj_coef = vec![50., 120.]; // c^T
    let cols: Vec<Col> = obj_coef.clone()
        .into_iter()
        .map(
            |x| {
                let col = m.add_col();
                if false {
                    m.set_integer(col);
                }
                col
            },
        )
        .collect(); // x 
    for i in 0..obj_coef.len(){
        m.set_obj_coeff(cols[i],obj_coef[i]); // setting up the objective function 
    }
    let A = CsMat::new((3, 2),
                       vec![0,2,4,6],
                       vec![0, 1, 0, 1, 0,1],
                       vec![100.,200.,10.,30.,1.,1.]); // matrix A    
    let upper = vec![10000.,1200.,110.]; // b
    // Ax <= b
    for i in 0..(A.rows()){
        let row = m.add_row();
        m.set_row_upper(row,upper[i]);
        for j in A.proper_indptr()[i]..(A.proper_indptr()[i+1]){
            m.set_weight(row,cols[A.indices()[j]] , A.data()[j]);
        }
    } 
    // Set objective sense.
    m.set_obj_sense(Sense::Maximize);

    // Solve the problem. Returns the solution
    let sol = m.solve();
 
    println!("{:?}", sol.col(cols[0]));
    println!("{:?}", sol.col(cols[1]));

    let dim = 1;
    let maxdis = 2;
 
    let ringmetadata = exhact::matrix::RingMetadata{
     ringspec: RingSpec::Rational,
     identity_additive: Ratio::new(0, 1),
     identity_multiplicative: Ratio::new(1, 1)
    };

    let dismat = vec![  vec![0,  1,  2,  1],
                        vec![1,  0,  1,  2],
                        vec![2,  1,  0,  1],
                        vec![1,  2,  1,  0]  ];
    let chx = exhact::clique::CliqueComplex {
        // the distance/dissimilarity matrix
        dissimilarity_matrix: dismat, 
        // threshold to stop the filtration
        dissimilarity_value_max: maxdis, 
        // sets "safeguards" on dimension; we'll get warnings if we try to 
        // get boundary matrices in dimension higher than dim+1
        safe_homology_degrees_to_build_boundaries: (1..dim+1).collect(), 
        // set the default major dimension (for sparse matrices) to be row
        major_dimension: MajorDimension::Row, 
        // indicates we want Z/3Z coefficients
        ringmetadata: ringmetadata, 
        // don't worry about this
        simplex_count: Vec::new() 
    };

    let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);


    for i in 1..(1+1){
        //record as simplices
        let simplex_bar = factored_complex.simplex_barcode(i);
        let birth = &simplex_bar[0].0;
        let death = &simplex_bar[0].1;
        // loop over Sn
        let Fn = chx.keys_unordered_itr(i).filter(|s| s <= &death );
        // loop over S_{n+1}
        let Fn1 = chx.keys_unordered_itr(i+1).filter(|s| s <= &death && s>=&birth);
        // build oracle for the entire boundary matrix
        let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row,
                              exhact::chx::ChxTransformKind::Boundary);
        // create hashmaps to store keys to indices 
        let mut maj_2_index = HashMap::new();       
        let mut min_2_index = HashMap::new();   
        // initialize indices to be 0   
        let maj_index = 0;
        let min_index = 0;
        // create sparse matrix
        for edge in Fn{
            if !maj_2_index.contains_key(&edge) {
                maj_2_index.insert(edge.clone(), maj_index); 
                let maj_index = maj_index + 1;
            }
            let minor_fields = D.maj_itr(&edge);
            for minor_field in minor_fields{
                let tri = minor_field.0;
                let data = minor_field.1;
                if !min_2_index.contains_key(&tri) {
                    min_2_index.insert(tri.clone(), min_index); 
                    let min_index = min_index + 1;
                }
                println!("c");
                println!("{:?}", tri);
                println!("{:?}", data);
            }
        }
    }
}