//! Linear optimization problem.

use ndarray::ArrayView1;

use crate::matrix::coo::CooMat;
use crate::problem::minlp::ProblemMinlp;
use crate::problem::milp::ProblemMilp;
use crate::problem::nlp::ProblemNlp;

/// Linear optimization problem (Lp) of the form
/// ```ignore
/// minimize   c^T*x       
/// subject to a*x = b     : lambda
///            l <= x <= u : pi and mu.
/// ```                     
pub struct ProblemLp {
    base_milp: ProblemMilp,
    base_nlp: ProblemNlp,
}

impl ProblemLp {

    /// Creates a new linear optimization problem (Lp).
    pub fn new(c: Vec<f64>,
               a: CooMat<f64>,
               b: Vec<f64>,  
               l: Vec<f64>,
               u: Vec<f64>,
               x0: Option<Vec<f64>>) -> Self {
        let nx = c.len();
        let base_milp = ProblemMilp::new(c.clone(), 
                                         a.clone(), 
                                         b.clone(), 
                                         l.clone(), 
                                         u.clone(), 
                                         vec![false; nx], 
                                         x0.clone());
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
        let base_nlp = ProblemNlp::new(CooMat::from_nnz((nx, nx), 0),
                                       a.clone(), 
                                       b.clone(),
                                       CooMat::from_nnz((0, nx), 0),
                                       Vec::new(), 
                                       l.clone(), 
                                       u.clone(), 
                                       x0.clone(),
                                       eval_fn);
        Self {
            base_milp: base_milp,
            base_nlp: base_nlp,
        }
    }

    /// Initial point.
    pub fn x0(&self) -> Option<&[f64]> { self.base_milp.x0() }

    /// Objective function gradient.
    pub fn c(&self) -> &[f64] { self.base_milp.c() }
    
    /// Jacobian matrix of linear equality constraints.
    pub fn a(&self) -> &CooMat<f64> { self.base_milp.a() } 
    
    /// Right-hand-side vector of linear equality constraints.
    pub fn b(&self) -> &[f64] { self.base_milp.b() }
    
    /// Vector of optimization variable lower limits.
    pub fn l(&self) -> &[f64] { self.base_milp.l() }
    
    /// Vector of optimization variable upper limits.
    pub fn u(&self) -> &[f64] { self.base_milp.u() }
    
    /// Number of optimization variables.
    pub fn nx(&self) -> usize { self.c().len() }
    
    /// Number of linear equality cosntraints.
    pub fn na(&self) -> usize { self.b().len() }
    
    /// Returns a mutable reference to the problem cast as a Milp.
    pub fn as_mut_milp(&mut self) -> &mut ProblemMilp { &mut self.base_milp }
    
    /// Returns a mutable reference to the problem cast as a Minlp.
    pub fn as_mut_minlp(&mut self) -> &mut ProblemMinlp { self.base_milp.as_mut_minlp() } 
    
    /// Returns a mutable reference to the problem cast as an Nlp
    pub fn as_mut_nlp(&mut self) -> &mut ProblemNlp { &mut self.base_nlp }      
}