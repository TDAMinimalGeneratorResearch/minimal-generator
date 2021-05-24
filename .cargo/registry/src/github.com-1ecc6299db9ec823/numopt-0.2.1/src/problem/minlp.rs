//! Mixed-integer nonlinear optimization problem.

use crate::matrix::coo::CooMat; 
use crate::problem::base::ProblemEval;

/// Mixed-integer nonlinear optimization problem (Minlp) of the form
/// ```ignore
/// minimize   phi(x)
/// subject to a*x = b
///            f(x) = 0
///            l <= x <= u
///            p*x in integers,
/// ```
/// where p*x gives a subvector of x.
pub struct ProblemMinlp
{
    x0: Option<Vec<f64>>,

    phi: f64,
    gphi: Vec<f64>,
    hphi: CooMat<f64>,  // lower triangular
    
    a: CooMat<f64>,
    b: Vec<f64>,
    
    f: Vec<f64>,
    j: CooMat<f64>,
    h: Vec<CooMat<f64>>,
    hcomb: CooMat<f64>, // lower triangular
    
    l: Vec<f64>,
    u: Vec<f64>,
    
    p: Vec<bool>,
    
    eval_fn: ProblemEval,
}

impl ProblemMinlp {

    /// Creates new mixed-integer nonlinear optimization 
    /// problem (Minlp).
    pub fn new(hphi: CooMat<f64>, 
               a: CooMat<f64>, 
               b: Vec<f64>,
               j: CooMat<f64>,
               h: Vec<CooMat<f64>>,  
               l: Vec<f64>, 
               u: Vec<f64>, 
               p: Vec<bool>,
               x0: Option<Vec<f64>>,
               eval_fn: ProblemEval) -> Self {

        let nx = a.cols();
        let na = a.rows();
        let nf = j.rows();

        assert_eq!(hphi.cols(), nx);
        assert_eq!(hphi.rows(), nx);

        assert_eq!(a.cols(), nx);
        assert_eq!(a.rows(), na);
        assert_eq!(b.len(), na);

        assert_eq!(j.cols(), nx);
        assert_eq!(j.rows(), nf);
        assert_eq!(h.len(), nf);
        for hh in h.iter() {
            assert_eq!(hh.rows(), nx);
            assert_eq!(hh.cols(), nx);
        }

        assert_eq!(l.len(), nx);
        assert_eq!(u.len(), nx);

        assert_eq!(p.len(), nx);

        let mut k: usize = 0;
        let hcomb_nnz = h.iter().map(|h| h.nnz()).sum();
        let mut hcomb = CooMat::from_nnz((nx, nx), hcomb_nnz);
        for hh in h.iter() {
            for (row, col, _val) in hh.iter() {
                hcomb.set_row_ind(k, *row);
                hcomb.set_col_ind(k, *col);
                k += 1;
            }
        }
        
        Self {
            x0: x0,
            phi: 0.,
            gphi: vec![0.;nx],
            hphi: hphi,
            a: a,
            b: b,
            f: vec![0.;nf],
            j: j,
            h: h,
            hcomb: hcomb,
            l: l,
            u: u,
            p: p,
            eval_fn: eval_fn
        }
    }

    /// Initial point.
    pub fn x0(&self) -> Option<&[f64]> { 
        match &self.x0 { 
            Some(xx) => Some(&xx),
            None => None
        }
    }

    /// Objective function value.
    pub fn phi(&self) -> f64 { self.phi }

    /// Objective function gradient value.
    pub fn gphi(&self) -> &[f64] { &self.gphi }

    /// Objective function Hessian value (lower triangular part)
    pub fn hphi(&self) -> &CooMat<f64> { &self.hphi }

    /// Jacobian matrix of linear equality constraints.
    pub fn a(&self) -> &CooMat<f64> { &self.a } 

    /// Right-hand-side vector of linear equality constraints.
    pub fn b(&self) -> &[f64] { &self.b }

    /// Nonlinear equality constraint function value.
    pub fn f(&self) -> &[f64] { &self.f }

    /// Nonlinear equality constraint function Jacobian value.
    pub fn j(&self) -> &CooMat<f64> { &self.j } 

    /// Vector of nonlinear equality constraint function Hessian values
    /// (lower triangular parts).
    pub fn h(&self) -> &Vec<CooMat<f64>> { &self.h } 

    /// Linear combination of nonlinear equality constraint function Hessian values
    /// (lower triangular parts).
    pub fn hcomb(&self) -> &CooMat<f64> { &self.hcomb }

    /// Vector of optimization variable lower limits.
    pub fn l(&self) -> &[f64] { &self.l }

    /// Vector of optimization variable upper limits.
    pub fn u(&self) -> &[f64] { &self.u }

    /// Vector of boolean values indicating optimization variables that are constrained
    /// to be integers.
    pub fn p(&self) -> &[bool] { &self.p } 
    
    /// Number of optimization variables.
    pub fn nx(&self) -> usize { self.a().cols() }

    /// Number of linear equality constraints.
    pub fn na(&self) -> usize { self.b().len() }

    /// Number of nonlinear equality constraints.
    pub fn nf(&self) -> usize { self.f().len() }
    
    /// Function that evaluates objective function, nonlinear equality constraint
    /// functions, and their derivatives for a given vector of optimization variable values.
    pub fn evaluate(&mut self, x: &[f64]) -> () {
        (self.eval_fn)(&mut self.phi, 
                       &mut self.gphi,
                       &mut self.hphi,
                       &mut self.f,
                       &mut self.j,
                       &mut self.h,
                       x)
    }

    /// Function that forms a linear combination of nonlinear equality constraint
    /// function Hessians.
    pub fn combine_h(&mut self, nu: &[f64]) -> () {
        assert_eq!(self.nf(), nu.len());
        let mut k: usize = 0;
        let data = self.hcomb.data_mut();
        for (h, nuval) in self.h.iter().zip(nu.iter()) {
            for val in h.data().iter() {
                data[k] = (*nuval)*(*val);
                k += 1;
            }
        }    
    }
}
