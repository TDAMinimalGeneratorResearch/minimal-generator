//! Base optimization solver types and trait.

use std::fmt;
use simple_error::SimpleError;
use std::collections::HashMap;

use crate::problem::base::{Problem, ProblemSol};

/// Optimization solver status.
#[derive(Debug, PartialEq)]
pub enum SolverStatus {

    /// Optimization solver successfully solved problem 
    /// to the specified accuracy.
    Solved,

    /// Optimization solver determined that problem
    /// is infeasible.
    Infeasible,

    /// Optimization solver has unknown status.
    Unknown,

    /// Optimization problem terminated with error.
    Error,
}

/// Optimization solver parameter.
#[derive(Clone)]
pub enum SolverParam {

    /// Integer solver parameter.
    IntParam(i32),

    /// Float solver parameter.
    FloatParam(f64),

    /// String solver parameter. 
    StrParam(String),
}

/// A trait for optimization solvers.
pub trait Solver {

    /// Gets optimization solver parameter value.
    fn get_param(&self, name: &str) -> Option<&SolverParam> { self.get_params().get(name) }

    /// Gets reference of optimization solver parameters.
    fn get_params(&self) -> &HashMap<String, SolverParam>;

    /// Gets mutable reference of optimization solver parameters.
    fn get_params_mut(&mut self) -> &mut HashMap<String, SolverParam>;
    
    /// Sets optimization solver parameter.
    fn set_param(&mut self, name: &str, value: SolverParam) -> Result<(), SimpleError> { 
       
        let v = match self.get_params_mut().get_mut(name) {
            Some(x) => x,
            None => return Err(SimpleError::new("unknown parameter"))
        };

        *v = match ((*v).clone(), value) {
            (SolverParam::IntParam(_x), SolverParam::IntParam(y)) => { 
                SolverParam::IntParam(y) 
            },
            (SolverParam::FloatParam(_x), SolverParam::FloatParam(y)) => { 
                SolverParam::FloatParam(y)
            },
            (SolverParam::StrParam(_x), SolverParam::StrParam(y)) => { 
                SolverParam::StrParam(y)
            }, 
            _ => return Err(SimpleError::new("invalid parameter type"))
        };    

        Ok(())
    }

    /// Solves optimization problem.
    fn solve(&self, problem: &mut Problem) -> Result<(SolverStatus, ProblemSol), SimpleError>;
}

impl Eq for SolverStatus {}

impl fmt::Display for SolverStatus {

    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolverStatus::Error => write!(f, "error"),
            SolverStatus::Unknown => write!(f, "unknown"),
            SolverStatus::Solved => write!(f, "solved"),
            SolverStatus::Infeasible => write!(f, "infeasible")
        }
    }
}

