//! Cbc solver interface.

use std::fs::File;
use std::ffi::OsStr;
use tempfile::Builder;
use std::io::prelude::*;
use std::fs::remove_file;
use std::process::Command;
use simple_error::SimpleError;
use std::io::{self, BufReader};
use std::collections::HashMap;

use crate::solver::base::{Solver, 
                          SolverParam,
                          SolverStatus};
use crate::problem::base::{Problem,
                           ProblemSol};
use crate::problem::milp::{ProblemMilp,
                           ProblemMilpIO};

/// Interface to the optimization solver Cbc from COIN-OR 
/// that utilzes the command-line tool "cbc". 
///
/// The command-line tool "cbc" needs to be on the system path.
/// 
/// It can solve problems of type [ProblemLp](../../problem/lp/struct.ProblemLp.html) 
/// and [ProblemMilp](../../problem/milp/struct.ProblemMilp.html).
pub struct SolverCbcCmd {
    parameters: HashMap<String, SolverParam>,
}

impl SolverCbcCmd {

    // Creates solver instance.
    pub fn new() -> Self { 

        let mut parameters: HashMap<String, SolverParam> = HashMap::new();
        parameters.insert("logLevel".to_string(), SolverParam::IntParam(5));

        Self {
            parameters: parameters,
        } 
    }

    /// Reads cbc solver solution file.
    pub fn read_sol_file(fname: &str, 
                         p: &ProblemMilp, 
                         cbc: bool) -> io::Result<(SolverStatus, ProblemSol)> {
        println!("here11");
        println!("here11");
        println!("here11");
        let mut name: String;
        let mut dtype: String;
        let mut index: usize;
        let mut value: f64;
        let mut mul: f64;
        let mut status = SolverStatus::Error;
        let mut solution = ProblemSol::new(p.nx(), p.na(), 0);
        let f = match File::open(fname) {
            Ok(ff) => ff,
            Err(_e) => return Ok((status, solution))
        };
        let mut r = BufReader::new(f);
        let mut line = String::new();
        let e = io::Error::new(io::ErrorKind::Other, "bad solution file");
        println!("here11");
        // Status
        r.read_line(&mut line)?;
        match line.split_ascii_whitespace().next() {
            Some(s) => {
                if !cbc && s == "optimal" {
                    status = SolverStatus::Solved;
                }
                else if !cbc && s == "infeasible" {
                    status = SolverStatus::Infeasible;
                }
                else if cbc && s == "Optimal" {
                    status = SolverStatus::Solved;
                }
                else if cbc && s == "Infeasible" {
                    status = SolverStatus::Infeasible;
                }
            },
            None => {
                status = SolverStatus::Error;
            }
        }
        println!("here2");
        // Not solved
        if status != SolverStatus::Solved {
            return Ok((status, solution))
        }

        // Objective value
        if !cbc {
            r.read_line(&mut line)?;
        }

        // Results
        for l in r.lines() {
            line = l?;
            let mut iter = line.split_ascii_whitespace();
            iter.next();
            name = match iter.next() {
                Some(s) => s.to_string(),
                None => return Err(e)
            };
            value = match iter.next() {
                Some(s) => match s.parse() { Ok(f) => f, Err(_e) => return Err(e) },
                None => return Err(e)
            };
            mul = match iter.next() {
                Some(s) => match s.parse() { Ok(f) => f, Err(_e) => return Err(e) },
                None => return Err(e)
            };
            let mut name_iter = name.split('_');
            dtype = match name_iter.next() {
                Some(s) => s.to_string(),
                None => return Err(e)
            };
            index = match name_iter.next() {
                Some(s) => match s.parse() { Ok(n) => n, Err(_e) => return Err(e) },
                None => return Err(e)
            };

            // Variable
             if dtype == "x" {
                solution.x[index] = value;
                if mul > 0. {
                    solution.pi[index] = mul;
                }
                else {
                    solution.mu[index] = -mul;
                }
            }

            // Constraint
            else if dtype == "c" {
                solution.lam[index] = mul;
            }
            else {
                return Err(e);
            }
        }

        Ok((status, solution))
    }
}

impl Solver for SolverCbcCmd {

    fn get_params(&self) -> &HashMap<String, SolverParam> { &self.parameters }
    fn get_params_mut(&mut self) -> &mut HashMap<String, SolverParam> { &mut self.parameters }
    
    fn solve(&self, problem: &mut Problem) -> Result<(SolverStatus, ProblemSol), SimpleError> {

        // Get problem
        let p  = match problem {
            Problem::Milp(x) => x,
            Problem::Lp(x) => x.as_mut_milp(),
            _ => return Err(SimpleError::new("problem type not supported"))
        };

        // Input filename
        let input_file = Builder::new()
            .prefix("cbc")
            .suffix(".lp")
            .tempfile();
        let input_filename = match input_file {
            Ok(f) => f.path().file_name().and_then(OsStr::to_str).unwrap().to_string(),
            Err(_e) => return Err(SimpleError::new("failed to create input filename")),
        };

        // Output filename
        let output_file = Builder::new()
            .prefix("cbc")
            .suffix(".sol")
            .tempfile();
        let output_filename = match output_file {
            Ok(f) => f.path().file_name().and_then(OsStr::to_str).unwrap().to_string(),
            Err(_e) => return Err(SimpleError::new("failed to create output filename")),
        };

        // Write input file
        match p.write_to_lp_file(&input_filename) {
            Ok(()) => (),
            Err(_e) => {
                remove_file(&input_filename).ok();
                remove_file(&output_filename).ok();
                return Err(SimpleError::new("failed to write lp file"));
            }
        };

        // Parameters
        let log_level = match self.get_param("logLevel") {
            Some(SolverParam::IntParam(i)) => i,
            _ => return Err(SimpleError::new("unable to get parameter logLevel"))
        };

        // Call Cbc command
        match Command::new("cbc")
                      .args(&[&input_filename,
                              "logLevel",
                              format!("{}", log_level).as_ref(),
                              "printingOptions",
                              "all", 
                              "solve", 
                              "solution",
                              &output_filename])
                      .spawn()
                      .and_then(|mut cmd| cmd.wait()) {
            Ok(_s) => (),
            Err(_e) => {
                remove_file(&input_filename).ok();
                remove_file(&output_filename).ok();
                return Err(SimpleError::new("failed executing cbc command"));
            }
        }
        
        // Clean up input file
        remove_file(&input_filename).ok();

        // Read output file
        let (status, solution) = match Self::read_sol_file(&output_filename, &p, true) {
            Ok((s, sol)) => (s, sol),
            Err(_e) => {
                remove_file(&output_filename).ok();
                return Err(SimpleError::new("failed to read cbc solution file"))
            }
        };

        // Clean up output file
        remove_file(&output_filename).ok();
        println!("here1");
        // All good
        Ok((status, solution))
    }
}

#[cfg(test)]
mod tests {

    use serial_test::serial;

    use crate::matrix::coo::CooMat;
    use crate::problem::base::Problem;
    use crate::problem::lp::ProblemLp;
    use crate::problem::milp::ProblemMilp;
    use crate::solver::base::{Solver, SolverParam, SolverStatus};
    use crate::solver::cbc_cmd::SolverCbcCmd;
    use crate::assert_vec_approx_eq;

    #[test]
    #[serial]
    fn cbc_solve_milp() {

        // Sample problem 
        // min        -x0 - x1 
        // subject to -2*x0 +  2*x1 + x2 == 1
        //            -8*x0 + 10*x1 + x3 ==  13
        //            x2 <= 0
        //            x3 >= 0
        //            x0 integer
        //            x1 integer

        let mut p = Problem::Milp(ProblemMilp::new(
            vec![-1.,-1., 0., 0.],
            CooMat::new(
                (2, 4),
                vec![0,0,0,1,1,1],
                vec![0,1,2,0,1,3],
                vec![-2.,2.,1.,-8.,10.,1.]),
            vec![1.,13.],
            vec![-1e8,-1e8,-1e8,0.],
            vec![1e8,1e8,0.,1e8],
            vec![true, true, false, false],
            None,
        ));

        let mut s = SolverCbcCmd::new();
        s.set_param("logLevel", SolverParam::IntParam(0)).unwrap();

        let (status, solution) = s.solve(&mut p).unwrap();

        assert_eq!(status, SolverStatus::Solved);
        assert_vec_approx_eq!(solution.x, 
                              &vec![1., 2., -1., 1.0], 
                              epsilon=1e-8);
    }

    #[test]
    #[serial]
    fn cbc_solve_lp() {

        // Sample problem 
        // min        180*x0 + 160*x1 
        // subject to 6*x0 +   x1 + x2 == 12
        //            3*x0 +   x1 + x3 ==  8
        //            4*x0 + 6*x1 + x4 == 24
        //            0 <= x0 <= 5
        //            0 <= x1 <= 5
        //            x2 <= 0
        //            x3 <= 0
        //            x4 <= 0

        let mut p = Problem::Lp(ProblemLp::new(
            vec![180.,160., 0., 0., 0.],
            CooMat::new(
                (3, 5),
                vec![0,0,0,1,1,1,2,2,2],
                vec![0,1,2,0,1,3,0,1,4],
                vec![6.,1.,1.,3.,1.,1.,4.,6.,1.]),
            vec![12.,8.,24.],
            vec![0.,0.,-1e8,-1e8,-1e8],
            vec![5.,5.,0.,0.,0.],
            None,
        ));

        let mut s = SolverCbcCmd::new();
        s.set_param("logLevel", SolverParam::IntParam(0)).unwrap();
        let (status, solution) = s.solve(&mut p).unwrap();

        assert_eq!(status, SolverStatus::Solved);
        assert_vec_approx_eq!(solution.x, 
                              &vec![1.7142857, 2.8571429, -1.1428571, 0., 0.], 
                              epsilon=1e-8);
        assert_vec_approx_eq!(solution.lam, 
                              &vec![0., 31.428571, 21.428571], 
                              epsilon=1e-8);
        assert_vec_approx_eq!(solution.mu, 
                              &vec![1.4210855e-14, 0., 0., 3.1428571e+01, 2.1428571e+01], 
                              epsilon=1e-8);
        assert_vec_approx_eq!(solution.pi, 
                              &vec![0.;5], 
                              epsilon=1e-8);

    }
}