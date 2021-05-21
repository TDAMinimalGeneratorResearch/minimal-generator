// 1. simplex wise barcode [key, matched_key] 

// 2. 

// 2. get index of key, matched_key: birth_ind, death_ind
// 3. for columns born after birth filtration, and dies before death_ind ==> M 
// 4. v[death_ind] == 1, len(v) = ncols(M)
// 5. M * v[birth_ind] != 0
// 6. M * v[birth_ind + 1, :] == 0 
// simplex_iterator.drop( |x|  has_wrong_birth(x)).collect();

use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind};
use exhact::clique::Simplex;
use std;

#[test]
fn main() {
    let dim = 3;
    let maxdis = 2 ;
    let ringmetadata = exhact::matrix::RingMetadata{
    	ringspec: RingSpec::Modulus(3),
    	identity_additive: 0,
    	identity_multiplicative: 1,
    };
    let dismat = vec![  vec![0,  1,  2,  1],
                        vec![1,  0,  1,  2],
                        vec![2,  1,  0,  1],
                        vec![1,  2,  1,  0]  ];
        
    
    // ----------------------------------------------------------------------------------
    // Construct the corresponding filtered clique complex
    // ----------------------------------------------------------------------------------
    let chx = exhact::clique::CliqueComplex {
        dissimilarity_matrix: dismat, 
        dissimilarity_value_max: maxdis, 
        safe_homology_degrees_to_build_boundaries: (1..dim+1).collect(), 
        major_dimension: MajorDimension::Row, 
        ringmetadata: ringmetadata, 
        simplex_count: Vec::new() 
    };

    let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row,
                              exhact::chx::ChxTransformKind::Boundary);
    
    let simplex_d1 = exhact::clique::Simplex{
        filvalue: 10,
        vertices : vec![0,1]
    };
    let simplex_d0 = exhact::clique::Simplex{
        filvalue: 0,
        vertices : vec![0]
    };
 
    let major_field = D.maj_itr( &simplex_d0 );

    let major_field2 = D.maj_itr( &simplex_d0 ); // this re-creates the iterator
 
    let minor_field = D.min_itr( &simplex_d1 );
 
    let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);
 
    let basis_vec = factored_complex.get_matched_basis_vector(1, &simplex_d1); 
 
    let basis_vec_iter = basis_vec
                            .iter() // get iterator for the hashmap
                            .map( |x| ( x.0.clone(), x.1.clone() ) ); // replace ref's with clones
    // println!("hhaa");
    // for item in basis_vec.iter()  {
    //     println!("{:?}", item);
    // };
    // get all simplices 
    let simplices = chx.keys_unordered_itr(1);
    
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
    println!("barrrrr");
    let barcode = factored_complex.barcode(1);
    for item in barcode.iter()  {
        println!("{:?}", item);
    }
    
};