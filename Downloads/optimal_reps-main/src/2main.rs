use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind};
use exhact::clique::Simplex;
use num::rational::Ratio;
use std;
use coin_cbc::{raw::Status, Col, Model, Sense, Solution as CbcSolution};
use sprs::{CsMat, CsVec};



fn main() {
    
// a simple knapsack problem
    // Maximize  5x[1] + 3x[2] + 2x[3] + 7x[4] + 4x[5]
    // s.t.      2x[1] + 8x[2] + 4x[3] + 2x[4] + 5x[5] <= 10
    //           
    // All x binary

    // Maximize  50x[1] + 20x[2]
    // s.t.      100x[1] + 200x[2]<= 10000
    //            10x[1] + 30x[2]<= 1200
    //              x[1] + x[2] <= 110
    // All x binary

    // Formulate an infeasible problem and try to solve it
    // Create the problem.
    let mut m = Model::default();
    m.set_parameter("log", "0");

    let c = vec![-5.0, -3.0];
    let cols: Vec<Col> = c
        .into_iter()
        .map(
            |x| {
                let col = m.add_col();
                if true {
                    m.set_integer(col);
                }
                col
            },
        )
        .collect();

    
    let obj_coef = vec![50., 120.];
    for i in 0..obj_coef.len(){
        m.set_obj_coeff(cols[i],obj_coef[i]);
    }

    // let mut A: CooMat<f64> = CooMat::<f64>::new(
    //     (3 , 2),
    //     vec![0, 0, 1, 1, 2,2],
    //     vec![0, 1, 0, 1, 0,1],
    //     vec![100.,200.,10.,30.,1.,1.]
    // );

    let A = CsMat::new((3, 2),
                       vec![0,2,4,6],
                       vec![0, 1, 0, 1, 0,1],
                       vec![100.,200.,10.,30.,1.,1.]);


    let upper = vec![10000.,1200.,110.];
    for i in 0..(A.rows()){
        let row = m.add_row();
        m.set_row_upper(row,upper[i]);
        for j in A.proper_indptr()[i]..(A.proper_indptr()[i+1]){
            m.set_weight(row,cols[A.indices()[j]] , A.data()[j]);
        }
    }





    // Set objective sense.
    m.set_obj_sense(Sense::Minimize);

    // Solve the problem. Returns the solution
    let sol = m.solve();

    // Check the result. sol.raw() returns a shared reference to the
    // raw bindings, allowing to use all getters.
    // assert_eq!(Status::Finished, sol.raw().status());
    // assert_eq!(5400., sol.raw().obj_value());
    println!("{:?}", sol.col(cols[0]));
    println!("{:?}", sol.col(cols[1]));
    // Check for the solution.
    // assert_eq!(60., sol.col(cols[0]));
    // assert_eq!(20., sol.col(cols[1]));




    let dim = 1;
    let maxdis = 2;
    
    // ----------------------------------------------------------------------------------
    // Define an object to represent the ring Z/3Z
    // ----------------------------------------------------------------------------------
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


    for i in 0..(dim+1){

        //record as simplices
        let simplex_bar = factored_complex.simplex_barcode(1);
        let birth = &simplex_bar[0].0;
        let death = &simplex_bar[0].1;
        println!("{:?}", simplex_bar);
        // loop over Sn
        let Fn = chx.keys_unordered_itr(1).filter(|s| s <= &death );
        // loop over S_{n+1}
        let Fn1 = chx.keys_unordered_itr(1+1).filter(|s| s <= &death && s>=&birth);
        // build oracle for the entire boundary matrix
        let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row,
                              exhact::chx::ChxTransformKind::Boundary);
        for item in Fn{
            let entries = D.maj_itr(&item);
            for item2 in entries{
                println!("{:?}",item2);
            }
            
        }
    }





    // let mut x: Option<Vec<f64>> = Some(vec![1.0,0.0,0.0,1.0,1.0]);
    // let mut c = vec![-5.0, -3.0, -2.0, -7.0, -4.0];
    // let mut b = vec![-9.0];
    // let mut u = vec![0.0, 0.0, 0.0, 0.0, 0.0];
    // let mut l = vec![-2.0,2.0,-2.0,-2.0,-2.0];
    // let mut p = vec![false, false, false, false, false];
    // let mut problem: numopt::problem::milp::ProblemMilp =  numopt::problem::milp::ProblemMilp::new(
    //     c,
    //     A,
    //     b,
    //     l,
    //     u,
    //     p,
    //     x
    // ); 
    // let sol = numopt::solver::cbc_cmd::SolverCbcCmd::read_sol_file(
    //     "file", 
    //     &problem, 
    //     true);
    // for item in sol.iter(){
    //     println!("{:?}", item);
    // }


    // // ----------------------------------------------------------------------------------
    // // Set maximum threshold values for homology dimension and dissimilarity
    // // ----------------------------------------------------------------------------------
    // let dim = 3;
    // let maxdis = 2;
    
    
    // // ----------------------------------------------------------------------------------
    // // Define an object to represent the ring Z/3Z
    // // ----------------------------------------------------------------------------------
    // let ringmetadata = exhact::matrix::RingMetadata{
    //  ringspec: RingSpec::Rational,
    //  identity_additive: Ratio::new(0, 1),
    //  identity_multiplicative: Ratio::new(1, 1)
    // };
    

    // // ----------------------------------------------------------------------------------
    // // Build a "dissimilarity matrix" as a vector of vectors
    // // ----------------------------------------------------------------------------------
    // let dismat = vec![  vec![0,  1,  2,  1],
    //                     vec![1,  0,  1,  2],
    //                     vec![2,  1,  0,  1],
    //                     vec![1,  2,  1,  0]  ];
    
    
    // // ----------------------------------------------------------------------------------
    // // Construct the corresponding filtered clique complex
    // // ----------------------------------------------------------------------------------
    // let chx = exhact::clique::CliqueComplex {
    //     // the distance/dissimilarity matrix
    //     dissimilarity_matrix: dismat, 
    //     // threshold to stop the filtration
    //     dissimilarity_value_max: maxdis, 
    //     // sets "safeguards" on dimension; we'll get warnings if we try to 
    //     // get boundary matrices in dimension higher than dim+1
    //     safe_homology_degrees_to_build_boundaries: (1..dim+1).collect(), 
    //     // set the default major dimension (for sparse matrices) to be row
    //     major_dimension: MajorDimension::Row, 
    //     // indicates we want Z/3Z coefficients
    //     ringmetadata: ringmetadata, 
    //     // don't worry about this
    //     simplex_count: Vec::new() 
    // };
    

    // // ----------------------------------------------------------------------------------
    // // Get a (row-major) sparse matrix oracle for the boundary operator
    // // ----------------------------------------------------------------------------------
    // let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row,
    //                           exhact::chx::ChxTransformKind::Boundary);
    
   
    // // ----------------------------------------------------------------------------------
    // // Define a weighted 0-simplex and 1-simplex
    // // ----------------------------------------------------------------------------------
    // let simplex_d1 = exhact::clique::Simplex{
    //     filvalue: 1,
    //     vertices : vec![0, 1]
    // };
    // let simplex_d0 = exhact::clique::Simplex{
    //     filvalue: 0,
    //     vertices : vec![0]
    // };


    // // ----------------------------------------------------------------------------------
    // // Compute / check the filtration values
    // // ----------------------------------------------------------------------------------
    // std::assert_eq!(chx.key_2_filtration( &simplex_d0 ), 0);
    // std::assert_eq!(chx.key_2_filtration( &simplex_d1 ), 1);


    // // ----------------------------------------------------------------------------------
    // // Access a row of the boundary matrix 
    // //  - the boundary matrix is "row-major" so we call rows "major fields"
    // // ----------------------------------------------------------------------------------
    // // create an iterator that runs over the structural nonzero entries of the row
    // // each item returned by the iterator is a tuple of form (weighted_simplex, coefficient)
    // let major_field = D.maj_itr( &simplex_d0 );

    
    // // the following for-loop should print the following:
    // // ```text
    // // Structural nonzero entries of a row corresponding to a 0-simplex:
    // // (Simplex { filvalue: 1, vertices: [0, 1] }, -1)
    // // (Simplex { filvalue: 2, vertices: [0, 2] }, -1)
    // // (Simplex { filvalue: 1, vertices: [0, 3] }, -1)
    // // ```
    // println!("Structural nonzero entries of a row corresponding to a 0-simplex:");
    // for item in major_field  {
    //     println!("{:?}", item);
    // }
    
    // // check to ensure that the output is correct:
    // let mut correct_val : Vec< (Simplex<i64>, i16) >  = 
    //                   vec![ (Simplex{ filvalue: 1, vertices: vec![0, 1] }, -1),
    //                         (Simplex{ filvalue: 2, vertices: vec![0, 2] }, -1),
    //                         (Simplex{ filvalue: 1, vertices: vec![0, 3] }, -1) ];
    
    // let major_field2 = D.maj_itr( &simplex_d0 ); // this re-creates the iterator
    // std::assert_eq!( major_field2.eq(correct_val.iter().map( |x| x.clone() ) ), true);
    

    // // ----------------------------------------------------------------------------------
    // // Access a column of the boundary matrix 
    // //  - the boundary matrix is "row-major" so we call columns "column fields"
    // // ----------------------------------------------------------------------------------
    // // create an iterator that runs over the structural nonzero entries of the row
    // // each item returned by the iterator is a tuple of form (weighted_simplex, coefficient)
    // let minor_field = D.min_itr( &simplex_d1 );
    
    // // the following for-loop should print the following:
    // // ```text
    // // Structural nonzero entries of a column corresponding to a 1-simplex:
    // // (Simplex { filvalue: 0, vertices: [1] }, 1)
    // // (Simplex { filvalue: 0, vertices: [0] }, -1)
    // // ```
    // println!("Structural nonzero entries of a column corresponding to a 1-simplex:");
    // for item in minor_field  {
    //     println!("{:?}", item);
    // }

    // // check to ensure that the output is correct:
    // let mut correct_val : Vec< (Simplex<i64>, i16) >  = 
    //                   vec![  (Simplex{ filvalue: 0, vertices: vec![1] },  1),
    //                          (Simplex{ filvalue: 0, vertices: vec![0] }, -1)  ];
    
    // let minor_field2 = D.min_itr( &simplex_d1 ); // this re-creates the iterator
    // std::assert_eq!( minor_field2.eq(correct_val.iter().map( |x| x.clone() ) ), true);


    // // ----------------------------------------------------------------------------------
    // // Factor the complex 
    // // ----------------------------------------------------------------------------------
    // let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);


    // // ----------------------------------------------------------------------------------
    // // Read barcodes / check for correctness
    // // ----------------------------------------------------------------------------------
    // // predefine the set of correct solutions
    // let correct_barcodes = vec![    vec![ (0, 2), (0, 1), (0, 1), (0, 1) ], // dimension 0
    //                                 vec![ (1, 2) ],                         // dimension 1
    //                                 vec![],                                 // dimension 2
    //                                 vec![]                                  // dimension 3
    //                             ];
   
    // // confirm that correct solutions and actual solutions are the same
    // for i in 0..dim+1 {
    //     assert_eq!( correct_barcodes[i], factored_complex.barcode(i) );
    // }


    // // ----------------------------------------------------------------------------------
    // // Get the cycle representative associated with a (weighted) simplex
    // //    - this may corresopnd to a "length-0 bar", depending on choice
    // // ----------------------------------------------------------------------------------
    // // get the vector *represented as a hashmap*
    // let basis_vec = factored_complex.get_matched_basis_vector(1, &simplex_d1);

    // // this for-loop should print the following:
    // // ```text
    // // Basis vector corresponding to a 1-simplex: 
    // // (Simplex { filvalue: 1, vertices: [0, 1] }, 1) 
    // // ```
    // println!("Basis vector corresponding to a 1-simplex:");
    // for item in basis_vec.iter()  {
    //     println!("{:?}", item);
    // }

    // // check to ensure that the output is correct:
    //     // write a vector with correct values
    // let mut correct_val : Vec< (Simplex<i64>, i16) >  = 
    //                   vec![  (Simplex{ filvalue: 1, vertices: vec![0, 1] },  1)  ];
    
    //     // convert the hashmap to an iterator  
        
    //     // remark: if we had tried to create this iterator via a command
    //     // > factored_complex.get_matched_basis_vector(1, &simplex_d1).iter().map( etc. )
    //     // then we would have gotten an error; the reason is that by declaring the `basis_vec`
    //     // variable name we signaled rust that it should have a certain lifetime; see 
    //     // https://stackoverflow.com/questions/54056268/temporary-value-is-freed-at-the-end-of-this-statement/54056716#54056716
    //     // for further discussion
    // let basis_vec_iter = basis_vec
    //                         .iter() // get iterator for the hashmap
    //                         .map( |x| ( x.0.clone(), x.1.clone() ) ); // replace ref's with clones

    //     // check equality
    // std::assert_eq!( basis_vec_iter.eq( correct_val.iter().cloned() ) , true);


    // // ----------------------------------------------------------------------------------
    // // Get the birth and death filtraitons of a chain 
    // //    - here "chain" is formalized as a hashmap mapping keys to coefficients
    // // ----------------------------------------------------------------------------------
    // // this part is not implemented yet.
}