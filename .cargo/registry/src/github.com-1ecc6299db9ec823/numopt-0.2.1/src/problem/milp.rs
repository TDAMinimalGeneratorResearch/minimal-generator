//! Mixed-integer linear optimization problem.

use std::fs::File;
use std::io::{self, Write, BufWriter};
use ndarray::ArrayView1;

use crate::matrix::coo::CooMat;
use crate::problem::minlp::ProblemMinlp;

/// Mixed-integer linear optimization problem (Milp) of the form
/// ```ignore
/// minimize   c^T*x
/// subject to a*x = b
///            l <= x <= u
///            p*x in integers,
/// ``` 
/// where p*x gives a subvector of x.
pub struct ProblemMilp {
    c: Vec<f64>,
    base: ProblemMinlp,
}

/// A trait for reading and writing mixed-integer linear 
/// optimization problems (Milp).
pub trait ProblemMilpIO {

    /// Reads problem from LP file.
    fn read_from_lp_file(filename: &str) -> io::Result<ProblemMilp>;

    /// Writes problem to LP file.
    fn write_to_lp_file(&self, filename: &str) -> io::Result<()>;
}

impl ProblemMilp {

    /// Creates new mixed-integer linear optimization
    /// problem (Milp).
    pub fn new(c: Vec<f64>,
               a: CooMat<f64>,
               b: Vec<f64>,  
               l: Vec<f64>,
               u: Vec<f64>, 
               p: Vec<bool>,
               x0: Option<Vec<f64>>) -> Self {
        let cc = c.clone();
        println!("here!!");
        let eval_fn = Box::new(move | phi: &mut f64, 
                                      gphi: &mut Vec<f64>, 
                                      _hphi: &mut CooMat<f64>,
                                      _f: &mut Vec<f64>,
                                      _j: &mut CooMat<f64>,
                                      _h: &mut Vec<CooMat<f64>>,
                                      x: &[f64] | {
            *phi = ArrayView1::from(&c).dot(&ArrayView1::from(x));
            gphi.copy_from_slice(&c);
        });
        let nx = a.cols();
        let base = ProblemMinlp::new(CooMat::from_nnz((nx, nx), 0), // Hphi
                                     a, 
                                     b,
                                     CooMat::from_nnz((0, nx), 0),  // J
                                     Vec::new(), 
                                     l, 
                                     u, 
                                     p, 
                                     x0,
                                     eval_fn);
        Self {
            c: cc,
            base: base,
        }
    }

    /// Initial point.
    pub fn x0(&self) -> Option<&[f64]> { self.base.x0() }

    /// Objective function gradient.
    pub fn c(&self) -> &[f64] { &self.c }

    /// Jacobian matrix of linear equality constraints.
    pub fn a(&self) -> &CooMat<f64> { &self.base.a() } 

    /// Right-hand-side vector of linear equality constraints.
    pub fn b(&self) -> &[f64] { &self.base.b() }

    /// Vector of optimization variable lower limits.
    pub fn l(&self) -> &[f64] { &self.base.l() }

    /// Vector of optimization variable upper limits.
    pub fn u(&self) -> &[f64] { &self.base.u() }

    /// Vector of boolean values indicating optimization variables that are constrained
    /// to be integers.
    pub fn p(&self) -> &[bool] { self.base.p() }

    /// Number of optimization variables.
    pub fn nx(&self) -> usize { self.c().len() }

    /// Number of linear equality cosntraints.
    pub fn na(&self) -> usize { self.b().len() }

    /// Returns a mutable reference to the problem cast as a Minlp.
    pub fn as_mut_minlp(&mut self) -> &mut ProblemMinlp { &mut self.base }
    
}

impl ProblemMilpIO for ProblemMilp {
    
    fn read_from_lp_file(_filename: &str) -> io::Result<ProblemMilp> {

        Err(io::Error::new(io::ErrorKind::Other, "not implemented"))
    }

    fn write_to_lp_file(&self, filename: &str) -> io::Result<()> {

        let mut pre: char;
        let mut j: usize;
        let mut d: f64;
        let mut b: f64;

        let f = File::create(filename)?;

        let mut w = BufWriter::new(f);

        // Objective
        w.write("Minimize\n".as_bytes())?;
        w.write(" obj:\n".as_bytes())?;
        for (i, c) in self.c().iter().enumerate() {
            if *c > 0. {
                pre = '+';
            }
            else if *c < 0. {
                pre = '-';
            }
            else {
                continue;
            }
            if c.abs() == 1. {
                w.write(format!("     {} x_{}\n", pre, i).as_bytes())?;
            }
            else {
                w.write(format!("     {} {:.10e} x_{}\n", 
                                pre, 
                                c.abs(), i).as_bytes())?;
            }
        }

        // Constraints
        w.write("Subject to\n".as_bytes())?;
        let mut a = self.a().to_csr();
        a.sum_duplicates();
        for i in 0..a.rows() {
            b = self.b()[i];
            w.write(format!("  c_{}:\n", i).as_bytes())?;
            for k in a.indptr()[i]..a.indptr()[i+1] {
                j = a.indices()[k];
                d = a.data()[k];
                if d > 0. {
                    pre = '+';
                }
                else if d < 0. {
                    pre = '-';
                }
                else {
                    continue;
                }
                if d.abs() == 1. {
                    w.write(format!("     {} x_{}\n", pre, j).as_bytes())?;
                }
                else {
                    w.write(format!("     {} {:.10e} x_{}\n", 
                                    pre, 
                                    d.abs(), 
                                    j).as_bytes())?;
                }
            }
            w.write(format!("     = {:.10e}\n", b).as_bytes())?;
        }

        // Bounds
        w.write("Bounds\n".as_bytes())?;
        for i in 0..self.nx() {
            w.write(format!(" {:.10e} <= x_{} <= {:.10e}\n",
                            self.l()[i],
                            i,
                            self.u()[i]).as_bytes())?;
        }

        // General
        w.write("General\n".as_bytes())?;
        for (i,f) in self.p().iter().enumerate() {
            if *f {
                w.write(format!(" x_{}\n", i).as_bytes())?;
            }
        }

        // End
        w.write("End\n".as_bytes())?;

        w.flush()?;

        Ok(())
    }
}