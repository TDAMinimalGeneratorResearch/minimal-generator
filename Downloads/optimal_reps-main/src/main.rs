// ================================================================
use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension, InvMod};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind};
use exhact::clique::Simplex;
use std;
use std::collections::BinaryHeap;
use std::cmp::Reverse;
use numopt::problem::lp::ProblemLp;
use numopt::matrix::coo::CooMat;
use good_lp::{constraint, default_solver, Solution, SolverModel, variables};
use std::error::Error;


// A   = [2 8 4 2 5]
// x   = [x x x x x]^T
// c^T = [5 3 2 7 4]
// -2x = 10
// minimize x
fn main() {
    let mut A: CooMat<f64> = CooMat::<f64>::new(
        (1 , 1 ),
        vec![0],
        vec![0],
        vec![-2.0]
    );
    
    let mut x: Option<Vec<f64>> = Some(vec![1.0]);
    let mut c = vec![1.0];
    let mut b = vec![10.0];
    let mut u = vec![0.0];
    let mut l = vec![-15.0];
    let mut p = vec![false ];
    let mut problem: numopt::problem::milp::ProblemMilp =  numopt::problem::milp::ProblemMilp::new(
        c,
        A,
        b,
        l,
        u,
        p,
        x
    ); 
    let sol = numopt::solver::cbc_cmd::SolverCbcCmd::read_sol_file(
        "file", 
        &problem, 
        true);
    for item in sol.iter(){
        println!("{:?}", item);
    }
    


    // let dim = 3;
    // let maxdis = 2.0 ;
    // let ringmetadata = exhact::matrix::RingMetadata{
    // 	ringspec: RingSpec::Modulus(3),
    // 	identity_additive: 0.0,
    // 	identity_multiplicative: 1.0,
    // };
    // let dismat = vec![  vec![0.0,  1.0,  2.0,  1.0],
    //                     vec![1.0,  0.0,  1.0,  2.0],
    //                     vec![2.0,  1.0,  0.0,  1.0],
    //                     vec![1.0,  2.0,  1.0,  0.0]  ];
          
    // let chx = exhact::clique::CliqueComplex {
    //     dissimilarity_matrix: dismat, dissimilarity_value_max: maxdis, safe_homology_degrees_to_build_boundaries: (1..dim+1).collect(), major_dimension: MajorDimension::Row, ringmetadata: ringmetadata, simplex_count: Vec::new() };

    // let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row,
    //                           exhact::chx::ChxTransformKind::Boundary);
 
    // // let major_field = D.maj_itr( &simplex_d0 );

    // // let major_field2 = D.maj_itr( &simplex_d0 ); // this re-creates the iterator
 
    // // let minor_field = D.min_itr( &simplex_d1 );
 
    // let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);
 
    // let basis_vec = factored_complex.get_matched_basis_vector(1, &simplex_d1); 
 
    // let basis_vec_iter = basis_vec
                            // .iter() // get iterator for the hashmap
                            // .map( |x| ( x.0.clone(), x.1.clone() ) ); // replace ref's with clones
    // println!("hhaa");
    // for item in basis_vec.iter()  {
    //     println!("{:?}", item);
    // };
    // get all simplices 
    
    
    // for key in chx.keys_unordered_itr(1) { // for each k simplex
    //     if factored_complex.dim_indexing[1].minkey_2_index.contains_key(&key) { continue; } // if this is a pivot column in d1, continue
    //     else if factored_complex.dim_indexing[2].majkey_2_index.contains_key(&key) {  // if this is a pivot row in d2
    //         let ind = factored_complex.dim_indexing[2].majkey_2_index[&key];          // d2 row index 
    //         let matched_key = &factored_complex.dim_indexing[2].index_2_minkey[ind];  // d2 column index 
    //         let diam1 = chx.key_2_filtration(&key); // d1 birth time 
    //         let diam2 = chx.key_2_filtration(matched_key); // d2 death time 
    //         if diam1 == diam2 { continue; } // if equal, not a feature
    //         let basis = factored_complex.get_matched_basis_vector(1, &key); 
    //         println!("ahhaha");
    //         for item in basis_vec.iter()  {
    //             println!("{:?}", item);
    //             println!("{:?}", diam1);
    //             println!("{:?}", diam2);
    //         }
    //     } else {
    //         let diam1 = chx.key_2_filtration(&key); // born 
    //         let diam2 = chx.max_filtration(); // never dies 
    //         if diam1 == diam2 { continue; } 
    //         let basis = factored_complex.get_matched_basis_vector(1, &key); 
    //         println!("hehehe");
    //         for item in basis_vec.iter()  {
    //             println!("{:?}", item);
    //         }
    //     }
    // }
    // println!("barrrrr");
    // let barcode = factored_complex.barcode(1);
    // for (birth, death) in barcode.iter()  { 
    //     // let birth_fil = chx.key_2_filtration( &birth );
    //     // let death_fil = chx.key_2_filtration( &death);
    //     println!("{:?}", birth);
    //     println!("{:?}", death);
        
    //     let tri_iterator = chx.keys_unordered_itr(2);
    //     let mut heap_tri = BinaryHeap::new();  
    //     let filtered_tri = tri_iterator.filter( |x| &x <= &death && &x >= &birth);
    //     for item in filtered_tri {
    //         heap_tri.push(Reverse(item.clone()))
    //     }

    //     let edge_iterator = chx.keys_unordered_itr(1);
    //     let mut heap_edge = BinaryHeap::new();  
    //     let filtered_edge = edge_iterator.filter( |x| &x <= &death && &x >= &birth);
    //     for item in filtered_edge {
    //         heap_edge.push(Reverse(item.clone()))
    //     }


    //     while let Some(Reverse(minkey)) = heap_edge.pop(){ // loop through columns 
    //         println!("{:?}", minkey);
    //     }
    // }
}