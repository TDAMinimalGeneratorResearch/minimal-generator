use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind};
use exhact::clique::Simplex;
use std;


fn main() {


    // ----------------------------------------------------------------------------------
    // Set maximum threshold values for homology dimension and dissimilarity
    // ----------------------------------------------------------------------------------
    let dim = 3;
    let maxdis = 2;
    
    
    // ----------------------------------------------------------------------------------
    // Define an object to represent the ring Z/3Z
    // ----------------------------------------------------------------------------------
    let ringmetadata = exhact::matrix::RingMetadata{
    	ringspec: RingSpec::Modulus(3),
    	identity_additive: 0,
    	identity_multiplicative: 1,
    };
    

    // ----------------------------------------------------------------------------------
    // Build a "dissimilarity matrix" as a vector of vectors
    // ----------------------------------------------------------------------------------
    let dismat = vec![  vec![0.0,  1.0,  2.0,  1.0],
                        vec![1.0,  0.0,  1.0,  2.0],
                        vec![2.0,  1.0,  0.0,  1.0],
                        vec![1.0,  2.0,  1.0,  0.0]  ];
    
    
    // ----------------------------------------------------------------------------------
    // Construct the corresponding filtered clique complex
    // ----------------------------------------------------------------------------------
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
    

    // ----------------------------------------------------------------------------------
    // Get a (row-major) sparse matrix oracle for the boundary operator
    // ----------------------------------------------------------------------------------
    let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row,
                              exhact::chx::ChxTransformKind::Boundary);
    
   
    // ----------------------------------------------------------------------------------
    // Define a weighted 0-simplex and 1-simplex
    // ----------------------------------------------------------------------------------
    let simplex_d1 = exhact::clique::Simplex{
        filvalue: 1,
        vertices : vec![0, 1]
    };
    let simplex_d0 = exhact::clique::Simplex{
        filvalue: 0,
        vertices : vec![0]
    };


    // ----------------------------------------------------------------------------------
    // Compute / check the filtration values
    // ----------------------------------------------------------------------------------
    std::assert_eq!(chx.key_2_filtration( &simplex_d0 ), 0);
    std::assert_eq!(chx.key_2_filtration( &simplex_d1 ), 1);


    // ----------------------------------------------------------------------------------
    // Access a row of the boundary matrix 
    //  - the boundary matrix is "row-major" so we call rows "major fields"
    // ----------------------------------------------------------------------------------
    // create an iterator that runs over the structural nonzero entries of the row
    // each item returned by the iterator is a tuple of form (weighted_simplex, coefficient)
    let major_field = D.maj_itr( &simplex_d0 );

    
    // the following for-loop should print the following:
    // ```text
    // Structural nonzero entries of a row corresponding to a 0-simplex:
    // (Simplex { filvalue: 1, vertices: [0, 1] }, -1)
    // (Simplex { filvalue: 2, vertices: [0, 2] }, -1)
    // (Simplex { filvalue: 1, vertices: [0, 3] }, -1)
    // ```
    // println!("Structural nonzero entries of a row corresponding to a 0-simplex:");
    // for item in major_field  {
    //     println!("{:?}", item);
    // }
    
    // check to ensure that the output is correct:
    let mut correct_val : Vec< (Simplex<i64>, i16) >  = 
                      vec![ (Simplex{ filvalue: 1, vertices: vec![0, 1] }, -1),
                            (Simplex{ filvalue: 2, vertices: vec![0, 2] }, -1),
                            (Simplex{ filvalue: 1, vertices: vec![0, 3] }, -1) ];
    
    let major_field2 = D.maj_itr( &simplex_d0 ); // this re-creates the iterator
    std::assert_eq!( major_field2.eq(correct_val.iter().map( |x| x.clone() ) ), true);
    

    // ----------------------------------------------------------------------------------
    // Access a column of the boundary matrix 
    //  - the boundary matrix is "row-major" so we call columns "column fields"
    // ----------------------------------------------------------------------------------
    // create an iterator that runs over the structural nonzero entries of the row
    // each item returned by the iterator is a tuple of form (weighted_simplex, coefficient)
    let minor_field = D.min_itr( &simplex_d1 );

    // the following for-loop should print the following:
    // ```text
    // Structural nonzero entries of a column corresponding to a 1-simplex:
    // (Simplex { filvalue: 0, vertices: [1] }, 1)
    // (Simplex { filvalue: 0, vertices: [0] }, -1)
    // ```
    // println!("Structural nonzero entries of a column corresponding to a 1-simplex:");
    // for item in minor_field  {
    //     println!("{:?}", item);
    // }
    
    // check to ensure that the output is correct:
    let mut correct_val : Vec< (Simplex<i64>, i16) >  = 
                      vec![  (Simplex{ filvalue: 0, vertices: vec![1] },  1),
                             (Simplex{ filvalue: 0, vertices: vec![0] }, -1)  ];
    
    let minor_field2 = D.min_itr( &simplex_d1 ); // this re-creates the iterator
    std::assert_eq!( minor_field2.eq(correct_val.iter().map( |x| x.clone() ) ), true);


    // ----------------------------------------------------------------------------------
    // Factor the complex 
    // ----------------------------------------------------------------------------------
    let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);


    // ----------------------------------------------------------------------------------
    // Read barcodes / check for correctness
    // ----------------------------------------------------------------------------------
    // predefine the set of correct solutions
    let correct_barcodes = vec![    vec![ (0, 2), (0, 1), (0, 1), (0, 1) ], // dimension 0
                                    vec![ (1, 2) ],                         // dimension 1
                                    vec![],                                 // dimension 2
                                    vec![]                                  // dimension 3
                                ];
   
    // confirm that correct solutions and actual solutions are the same
    for i in 0..dim+1 {
        assert_eq!( correct_barcodes[i], factored_complex.barcode(i) );
    }


    // ----------------------------------------------------------------------------------
    // Get the cycle representative associated with a (weighted) simplex
    //    - this may corresopnd to a "length-0 bar", depending on choice
    // ----------------------------------------------------------------------------------
    // get the vector *represented as a hashmap*
    let basis_vec = factored_complex.get_matched_basis_vector(1, &simplex_d1); 

    // this for-loop should print the following:
    // ```text
    // Basis vector corresponding to a 1-simplex: 
    // (Simplex { filvalue: 1, vertices: [0, 1] }, 1) 
    // ```
    // println!("Basis vector corresponding to a 1-simplex:");
    // for item in basis_vec.iter()  {
    //     println!("{:?}", item);
    // }

    // check to ensure that the output is correct:
        // write a vector with correct values
    let mut correct_val : Vec< (Simplex<i64>, i16) >  = 
                      vec![  (Simplex{ filvalue: 1, vertices: vec![0, 1] },  1)  ];
    
        // convert the hashmap to an iterator  
        
        // remark: if we had tried to create this iterator via a command
        // > factored_complex.get_matched_basis_vector(1, &simplex_d1).iter().map( etc. )
        // then we would have gotten an error; the reason is that by declaring the `basis_vec`
        // variable name we signaled rust that it should have a certain lifetime; see 
        // https://stackoverflow.com/questions/54056268/temporary-value-is-freed-at-the-end-of-this-statement/54056716#54056716
        // for further discussion
    let basis_vec_iter = basis_vec
                            .iter() // get iterator for the hashmap
                            .map( |x| ( x.0.clone(), x.1.clone() ) ); // replace ref's with clones

        // check equality
    std::assert_eq!( basis_vec_iter.eq( correct_val.iter().cloned() ) , true);


    // ----------------------------------------------------------------------------------
    // Get the birth and death filtraitons of a chain 
    //    - here "chain" is formalized as a hashmap mapping keys to coefficients
    // ----------------------------------------------------------------------------------
    // this part is not implemented yet.
}

