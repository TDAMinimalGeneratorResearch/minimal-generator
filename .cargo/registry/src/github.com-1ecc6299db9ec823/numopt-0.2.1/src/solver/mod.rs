//! Optimization solver interfaces.

pub mod base;
pub mod clp_cmd;
pub mod cbc_cmd;

#[cfg(feature = "ipopt")] 
pub mod ipopt;

pub use base::{Solver, SolverParam, SolverStatus};

#[cfg(feature = "ipopt")] 
pub use ipopt::SolverIpopt;

pub use clp_cmd::SolverClpCmd;
pub use cbc_cmd::SolverCbcCmd;