//! Numerical optimization problem abstractions, solver interfaces, and modeling tools.
//! 
//! ## Features
//! - Abstractions for Minlp, Nlp, Milp, and Lp optimization problems.
//! - Interfaces for COIN-OR optimization solvers Cbc, Clp, and Ipopt.
//! - Modeling tools with automatic sparse first- and second-order derivatives.

pub mod problem;
pub mod solver;
pub mod model;
pub mod matrix;
pub mod macros;