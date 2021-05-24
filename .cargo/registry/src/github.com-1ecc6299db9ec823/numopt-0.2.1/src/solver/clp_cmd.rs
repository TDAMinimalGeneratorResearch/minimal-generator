//! Clp solver interface.

use std::ffi::OsStr;
use tempfile::Builder;
use std::fs::remove_file;
use std::process::{Command, Stdio};
use simple_error::SimpleError;
use std::collections::HashMap;

use crate::solver::base::{Solver, 
                          SolverParam,
                          SolverStatus};
use crate::solver::cbc_cmd::SolverCbcCmd;
use crate::problem::base::{Problem,
                          ProblemSol};
use crate::problem::milp::ProblemMilpIO;

/// Interface to the optimization solver Clp from COIN-OR 
/// that utilzes the command-line tool "clp". 
/// 
/// The command-line tool "clp" needs to be on the system path.
/// 
/// It can solve problems of type [ProblemLp](../../problem/lp/struct.ProblemLp.html). 
pub struct SolverClpCmd {
    parameters: HashMap<String, SolverParam>,
}

impl SolverClpCmd {

    // Creates solver instance.
    pub fn new() -> Self { 

        let mut parameters: HashMap<String, SolverParam> = HashMap::new();
        parameters.insert("logLevel".to_string(), SolverParam::IntParam(5));

        Self {
            parameters: parameters,
        } 
    }
}

impl Solver for SolverClpCmd {

    fn get_params(&self) -> &HashMap<String, SolverParam> { &self.parameters }
    fn get_params_mut(&mut self) -> &mut HashMap<String, SolverParam> { &mut self.parameters }

    fn solve(&self, problem: &mut Problem) -> Result<(SolverStatus, ProblemSol), SimpleError> {

        // Get problem
        let p  = match problem {
            Problem::Lp(x) => x, //x.as_mut_milp(),
            _ => return Err(SimpleError::new("problem type not supported"))
        };

        // Input filename
        let input_file = Builder::new()
            .prefix("clp")
            .suffix(".lp")
            .tempfile();
        let input_filename = match input_file {
            Ok(f) => f.path().file_name().and_then(OsStr::to_str).unwrap().to_string(),
            Err(_e) => return Err(SimpleError::new("failed to create input filename")),
        };

        // Output filename
        let output_file = Builder::new()
            .prefix("clp")
            .suffix(".sol")
            .tempfile();
        let output_filename = match output_file {
            Ok(f) => f.path().file_name().and_then(OsStr::to_str).unwrap().to_string(),
            Err(_e) => return Err(SimpleError::new("failed to create output filename")),
        };

        // Write input file
        match p.as_mut_milp().write_to_lp_file(&input_filename) {
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

        // Call Clp command
        match Command::new("clp")
                      .stdout(if *log_level == 0 { Stdio::null() } else { Stdio::inherit() })
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
                return Err(SimpleError::new("failed executing clp command"));
            }
        }
        
        // Clean up input file
        remove_file(&input_filename).ok();

        // Read output file
        let (status, solution) = match SolverCbcCmd::read_sol_file(&output_filename, 
                                                                   p.as_mut_milp(), 
                                                                   false) {
            Ok((s, sol)) => (s, sol),
            Err(_e) => {
                remove_file(&output_filename).ok();
                return Err(SimpleError::new("failed to read clp solution file"))
            }
        };

        // Clean up output file
        remove_file(&output_filename).ok();
        
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
    use crate::solver::base::{Solver, SolverParam, SolverStatus};
    use crate::solver::clp_cmd::SolverClpCmd;
    use crate::assert_vec_approx_eq;

    #[test]
    #[serial]
    fn clp_solve_lp() {

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

        let mut s = SolverClpCmd::new();
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