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


// use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
// use exhact::cubical::{Cube, CubicalComplex};
// use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind, FactoredComplexBlockCsm, Indexing};
// use exhact::clique::Simplex;
// use std;
// use std::marker::PhantomData;
// use exhact::csm::CSM;
// use exhact::decomp::decomposition;


// fn main() {

//     // ----------------------------------------------------------------------------------
//     // Set maximum threshold values for homology dimension and dissimilarity
//     // ----------------------------------------------------------------------------------
//     let dim = 3;
//     let maxdis = 14;
    
    
//     // ----------------------------------------------------------------------------------
//     // Define an object to represent the ring Z/3Z
//     // ----------------------------------------------------------------------------------
//     let ringmetadata = exhact::matrix::RingMetadata{
//     	ringspec: RingSpec::Modulus(3),
//     	identity_additive: 0,
//     	identity_multiplicative: 1,
//     };
//     println!("haha");
//     // ----------------------------------------------------------------------------------
//     // Build a "dissimilarity matrix" as a vector of vectors
//     // ----------------------------------------------------------------------------------
//     let dismat = vec![  vec![0,  10,  14,  10],
//                         vec![10,  0,  10,  15],
//                         vec![14,  10,  0,  10],
//                         vec![10,  15,  10,  0]  ];
    
    
//     // ----------------------------------------------------------------------------------
//     // Construct the corresponding filtered clique complex
//     // ----------------------------------------------------------------------------------
//     let chx = exhact::clique::CliqueComplex {
//         // the distance/dissimilarity matrix
//         dissimilarity_matrix: dismat, 
//         // threshold to stop the filtration
//         dissimilarity_value_max: maxdis, 
//         // sets "safeguards" on dimension; we'll get warnings if we try to 
//         // get boundary matrices in dimension higher than dim+1
//         safe_homology_degrees_to_build_boundaries: (1..dim+1).collect(), 
//         // set the default major dimension (for sparse matrices) to be row
//         major_dimension: MajorDimension::Row, 
//         // indicates we want Z/3Z coefficients
//         ringmetadata: ringmetadata, 
//         // don't worry about this
//         simplex_count: Vec::new() 
//     };
//     let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row,
//                               exhact::chx::ChxTransformKind::Boundary);
    
//     let simplex_d1 = exhact::clique::Simplex{
//         filvalue: 1,
//         vertices : vec![3]
//     };
//     let simplex_d0 = exhact::clique::Simplex{
//         filvalue: 0,
//         vertices : vec![0]
//     };
 
//     let major_field = D.maj_itr( &simplex_d1 );
//     println!("Structural nonzero entries of a row corresponding to a 0-simplex:");
//     for item in major_field  {
//         println!("{:?}", item);
//     }
//     let major_field = D.maj_itr( &simplex_d0 );
//     println!("Structural nonzero entries of a row corresponding to a 0-simplex:");
//     for item in major_field  {
//         println!("{:?}", item);
//     }

//     std::assert_eq!(chx.key_2_filtration( &simplex_d0 ), 0);
//     std::assert_eq!(chx.key_2_filtration( &simplex_d1 ), 1);
 
//     // -------------------------------
//     let max_homology_degree = dim+1;
//     println!("here");
//     let original_complex = &chx;
//     let hmatrix = original_complex.get_smoracle(
//         MajorDimension::Row,
//         ChxTransformKind::Boundary
//     );

//     let mut blocks = FactoredComplexBlockCsm{
//         phantom: PhantomData,
//         original_complex: original_complex,
//         dim_rowoper: vec![CSM::new(MajorDimension::Row, hmatrix.ring().clone())],
//         dim_indexing: vec![Indexing::new()]
//     };

//     let mut maj_to_reduce = Vec::new();
//     for ii in 1..(max_homology_degree+1){
//         println!("dim{:?}", ii);
//         maj_to_reduce.clear();
//         let major_keys = original_complex.keys_ordered(ii-1);
//         for key in major_keys.iter() {
//             if !blocks.dim_indexing[ii-1].minkey_2_index.contains_key(key) {
//                 maj_to_reduce.push(key.clone());
//             }
//         }
//         let (rowoper, indexing) = decomposition(&hmatrix, &mut maj_to_reduce);
//         blocks.dim_rowoper.push(rowoper);
//         blocks.dim_indexing.push(indexing);
//     }


 
// }

//! In this example, we have multiple products,
//! and each consumes a set amount of fuel, and of time to produce.
//! The goal is to find, knowing the available fuel and time,
//! and the value of each product, how much we should produce of each.
//!
//! In this example, the number of resources is fixed (only fuel an time),
//! and the amount of products varies.
//! In the opposite case (a fixed number of products and an arbitrary number of resources),
//! the modelling is even simpler: you don't have to store any expression in your problem struct,
//! you can instantiate a SolverModel directly when creating your problem,
//! and then use SolverModel::with to add constraints dynamically.

// use good_lp::variable::ProblemVariables;
// use good_lp::{default_solver, variable, variables, Expression, Solution, SolverModel, Variable};

// struct Product {
//     // amount of fuel producing 1 unit takes
//     needed_fuel: f64,
//     // time it takes to produce 1 unit
//     needed_time: f64,
//     value: f64, // The amount of money we can sell an unit of the product for
// }

// struct ResourceAllocationProblem {
//     vars: ProblemVariables,
//     total_value: Expression,
//     consumed_fuel: Expression,
//     consumed_time: Expression,
//     available_fuel: f64,
//     available_time: f64,
// }

// impl ResourceAllocationProblem {
//     fn new(available_fuel: f64, available_time: f64) -> ResourceAllocationProblem {
//         ResourceAllocationProblem {
//             vars: variables!(),
//             available_fuel,
//             available_time,
//             consumed_fuel: 0.into(),
//             consumed_time: 0.into(),
//             total_value: 0.into(),
//         }
//     }

//     /// Add a new product to take into account in the optimization
//     fn add(&mut self, product: Product) -> Variable {
//         let amount_to_produce = self.vars.add(variable().min(0));
//         self.total_value += amount_to_produce * product.value;
//         self.consumed_fuel += amount_to_produce * product.needed_fuel;
//         self.consumed_time += amount_to_produce * product.needed_time;
//         amount_to_produce
//     }

//     fn best_product_quantities(self) -> impl Solution {
//         self.vars
//             .maximise(self.total_value)
//             .using(default_solver)
//             .with(self.consumed_fuel.leq(self.available_fuel))
//             .with(self.consumed_time.leq(self.available_time))
//             .solve()
//             .unwrap()
//     }
// }

// use float_eq::assert_float_eq;

// #[test]
// fn resource_allocation() {
//     let mut pb = ResourceAllocationProblem::new(5., 3.);
//     let steel = pb.add(Product {
//         needed_fuel: 1.,
//         needed_time: 1.,
//         value: 10.,
//     });
//     let stainless_steel = pb.add(Product {
//         needed_fuel: 2.,
//         needed_time: 1.,
//         value: 11.,
//     });

//     let solution = pb.best_product_quantities();

//     // The amount of steel we should produce
//     println!("solution value");
//     println!("{}", solution.value(steel));
//     println!("solution value done");
//     assert_float_eq!(1., solution.value(steel), abs <= 1e-10);
//     // The amount of stainless steel we should produce
//     assert_float_eq!(2., solution.value(stainless_steel), abs <= 1e-10);
// }
// fn main(){

// }

// #[test]
// fn using_a_vector() {
//     let products = vec![
//         Product {
//             needed_fuel: 1.,
//             needed_time: 1.,
//             value: 10.,
//         },
//         Product {
//             needed_fuel: 2.,
//             needed_time: 1.,
//             value: 11.,
//         },
//     ];

//     let mut pb = ResourceAllocationProblem::new(5., 3.);
//     let variables: Vec<_> = products.into_iter().map(|p| pb.add(p)).collect();
//     let solution = pb.best_product_quantities();
//     let product_quantities: Vec<_> = variables.iter().map(|&v| solution.value(v)).collect();
//     assert_float_eq!(1., product_quantities[0], abs <= 1e-10);
//     assert_float_eq!(2., product_quantities[1], abs <= 1e-10);
// }