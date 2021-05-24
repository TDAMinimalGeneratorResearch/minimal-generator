//! Ipopt solver interface.

#![cfg(feature = "ipopt")] 

use std::{ptr, slice};
use std::ffi::CString;
use ndarray::ArrayView1;
use simple_error::SimpleError;
use std::collections::HashMap;
use num_traits::cast::ToPrimitive;
use libc::{c_int, c_void, c_char, c_double};

#[path = "cipopt.rs"]
mod cipopt;

use crate::solver::base::{Solver, 
                          SolverParam,
                          SolverStatus};
use crate::problem::base::{Problem, ProblemSol};
use crate::problem::nlp::ProblemNlp;

/// Interface to the optimization solver Ipopt from COIN-OR
/// that links with the library "libipopt".
///  
/// The library needs to be on the linker path.
/// This solver is available when the feature "ipopt" is enabled.
/// 
/// It can solve problems of type [ProblemLp](../../problem/lp/struct.ProblemLp.html) 
/// and [ProblemNlp](../../problem/nlp/struct.ProblemNlp.html).
pub struct SolverIpopt {
    parameters: HashMap<String, SolverParam>,
}

impl SolverIpopt {

    // Creates solver instance.
    pub fn new() -> Self {
        
        let mut parameters: HashMap<String, SolverParam> = HashMap::new();
        parameters.insert("print_level".to_string(), SolverParam::IntParam(5));
        parameters.insert("sb".to_string(), SolverParam::StrParam(String::from("no")));
        
        Self {
            parameters: parameters,
        } 
    }
}

impl Solver for SolverIpopt {

    fn get_params(&self) -> &HashMap<String, SolverParam> { &self.parameters }
    fn get_params_mut(&mut self) -> &mut HashMap<String, SolverParam> { &mut self.parameters }
    
    fn solve(&self, problem: &mut Problem) -> Result<(SolverStatus, ProblemSol), SimpleError> {

        // Get problem
        let p  = match problem {
            Problem::Nlp(x) => x,
            Problem::Lp(x) => x.as_mut_nlp(),
            _ => return Err(SimpleError::new("problem type not supported"))
        };

        // Init
        let n: c_int = p.nx().to_i32().unwrap();
        let m: c_int = (p.na() + p.nf()).to_i32().unwrap();
        let glu: Vec<f64> = vec![0.; p.na()+p.nf()];
        let nnzj: c_int = (p.a().nnz() + p.j().nnz()).to_i32().unwrap();
        let nnzh: c_int = (p.hphi().nnz() + p.hcomb().nnz()).to_i32().unwrap();

        // Problem
        let cprob: cipopt::IpoptProblem = unsafe {
            cipopt::CreateIpoptProblem(n, 
                                       p.l().as_ptr(), 
                                       p.u().as_ptr(), 
                                       m, 
                                       glu.as_ptr(), 
                                       glu.as_ptr(), 
                                       nnzj, 
                                       nnzh, 
                                       0, 
                                       eval_f_cb, 
                                       eval_g_cb, 
                                       eval_grad_f_cb, 
                                       eval_jac_g_cb, 
                                       eval_h_cb)
        };
        if cprob.is_null() {
            return Err(SimpleError::new("failed to create ipopt problem"))
        }

        // Parameter print_level
        let print_level: &i32 = match self.get_param("print_level") {
            Some(SolverParam::IntParam(i)) => i,
            _ => return Err(SimpleError::new("unable to get parameter print_level"))
        };
        unsafe { add_int_option(cprob, "print_level", *print_level)? };

        // Parameter sb
        let sb: &str = match self.get_param("sb") {
            Some(SolverParam::StrParam(s)) => s.as_ref(),
            _ => return Err(SimpleError::new("unable to get parameter sb"))
        };
        unsafe { add_str_option(cprob, "sb", sb)? };
        
        // Buffers
        let mut x: Vec<f64> = vec![0.;p.nx()];
        let mut lamnu: Vec<f64> = vec![0.;p.na()+p.nf()];        
        let mut pi: Vec<f64> = vec![0.;p.nx()];
        let mut mu: Vec<f64> = vec![0.;p.nx()];

        // Initial point
        match p.x0() {
            Some(x0) => x.copy_from_slice(x0),
            None => ()
        };

        // Solve
        let cstatus : c_int = unsafe {
            cipopt::IpoptSolve(cprob, 
                               x.as_mut_ptr(), 
                               ptr::null_mut(), 
                               ptr::null_mut(), 
                               lamnu.as_mut_ptr(), 
                               pi.as_mut_ptr(), 
                               mu.as_mut_ptr(), 
                               p as *mut _ as *mut c_void)
        };

        // Set status and solution
        let mut status = SolverStatus::Error;
        if cstatus == 0 || cstatus == 1 {
            status = SolverStatus::Solved;
        }
        else if cstatus == 2 {
            status = SolverStatus::Infeasible;
        }  
        let mut solution = ProblemSol::new(p.nx(), p.na(), p.nf()); 
        solution.x.copy_from_slice(&x);
        for k in 0..p.na() {
            solution.lam[k] = -lamnu[k];
        }
        for k in p.na()..(p.na()+p.nf()) {
            solution.nu[k-p.na()] = -lamnu[k];
        }
        solution.pi.copy_from_slice(&pi);
        solution.mu.copy_from_slice(&mu);

        // Clean up
        unsafe {
            cipopt::FreeIpoptProblem(cprob);
        };

        Ok((status, solution))
    }
}

extern fn eval_f_cb(n: c_int, 
                    x: *const c_double, 
                    new_x: c_int, 
                    obj_value: *mut c_double, 
                    user_data: *mut c_void) -> c_int {
    unsafe {
        if x.is_null() || obj_value.is_null() || user_data.is_null(){
            return cipopt::FALSE;
        }
        let p: &mut ProblemNlp = &mut *(user_data as *mut ProblemNlp);
        match n.to_usize() {
            Some(nn) => { if nn != p.nx() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        if new_x == cipopt::TRUE {
            p.evaluate(slice::from_raw_parts(x, p.nx()));
        }
        *obj_value = p.phi();
    };
    cipopt::TRUE
}

extern fn eval_grad_f_cb(n: c_int, 
                         x: *const c_double, 
                         new_x: c_int, 
                         grad_f: *mut c_double, 
                         user_data: *mut c_void) -> c_int {
    unsafe {
        if x.is_null() || grad_f.is_null() || user_data.is_null() {
            return cipopt::FALSE;
        }
        let p: &mut ProblemNlp = &mut *(user_data as *mut ProblemNlp);
        match n.to_usize() {
            Some(nn) => { if nn != p.nx() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        if new_x == cipopt::TRUE {
            p.evaluate(slice::from_raw_parts(x, p.nx()));
        }
        ptr::copy(p.gphi().as_ptr(), grad_f, p.nx());
    };
    cipopt::TRUE
}

extern fn eval_g_cb(n: c_int, 
                    x: *const c_double, 
                    new_x: c_int, 
                    m: c_int,
                    g: *mut c_double, 
                    user_data: *mut c_void) -> c_int {
    unsafe {
        if x.is_null() || g.is_null() || user_data.is_null() {
            return cipopt::FALSE;
        }
        let p: &mut ProblemNlp = &mut *(user_data as *mut ProblemNlp);
        match n.to_usize() {
            Some(nn) => { if nn != p.nx() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        match m.to_usize() {
            Some(mm) => { if mm != p.na()+p.nf() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        let xx: &[f64] = slice::from_raw_parts(x, p.nx());
        if new_x == cipopt::TRUE {
            p.evaluate(&xx);
        }
        let ax = p.a()*xx.to_vec();
        let axmb = &ArrayView1::from(&ax)-&ArrayView1::from(p.b());
        ptr::copy(axmb.as_slice().unwrap().as_ptr(), g, p.na());
        ptr::copy(p.f().as_ptr(), g.add(p.na()), p.nf());
    };
    cipopt::TRUE
}

extern fn eval_jac_g_cb(n: c_int, 
                        x: *const c_double, 
                        new_x: c_int, 
                        m: c_int,
                        nele_jac: c_int,
                        irow: *mut c_int,
                        jcol: *mut c_int,
                        values: *mut c_double, 
                        user_data: *mut c_void) -> c_int {
    unsafe {
        if user_data.is_null() {
            return cipopt::FALSE;
        }
        let p: &mut ProblemNlp = &mut *(user_data as *mut ProblemNlp);
        match n.to_usize() {
            Some(nn) => { if nn != p.nx() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        match m.to_usize() {
            Some(mm) => { if mm != p.na()+p.nf() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        match nele_jac.to_usize() {
            Some(nnz) => { if nnz != p.a().nnz()+p.j().nnz() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        if values.is_null() {
            if irow.is_null() || jcol.is_null() {
                return cipopt::FALSE;
            }
            let mut k: usize = 0;
            for (row, col, _val) in p.a().iter() {
                *irow.add(k) = row.to_i32().unwrap();
                *jcol.add(k) = col.to_i32().unwrap();
                k += 1;
            }
            for (row, col, _val) in p.j().iter() {
                *irow.add(k) = row.to_i32().unwrap();
                *jcol.add(k) = col.to_i32().unwrap();
                k += 1;
            }
        }
        else {
            if new_x == cipopt::TRUE {
                if x.is_null() {
                    return cipopt::FALSE;
                }
                p.evaluate(slice::from_raw_parts(x, p.nx()));
            }
            let mut k: usize = 0;
            for (_row, _col, val) in p.a().iter() {
                *values.add(k) = *val; 
                k += 1;
            }
            for (_row, _col, val) in p.j().iter() {
                *values.add(k) = *val; 
                k += 1;
            }   
        }
    };
    cipopt::TRUE
}

extern fn eval_h_cb(n: c_int, 
                    x: *const c_double, 
                    new_x: c_int,
                    obj_factor: c_double, 
                    m: c_int,
                    lambda: *const c_double,
                    new_lambda: c_int,
                    nele_hess: c_int,
                    irow: *mut c_int,
                    jcol: *mut c_int,
                    values: *mut c_double, 
                    user_data: *mut c_void) -> c_int {
    unsafe {
        if user_data.is_null() {
            return cipopt::FALSE;
        }
        let p: &mut ProblemNlp = &mut *(user_data as *mut ProblemNlp);
        match n.to_usize() {
            Some(nn) => { if nn != p.nx() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        match m.to_usize() {
            Some(mm) => { if mm != p.na()+p.nf() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        match nele_hess.to_usize() {
            Some(nnz) => { if nnz != p.hphi().nnz()+p.hcomb().nnz() { return cipopt::FALSE; } },
            None => return cipopt::FALSE,
        }
        if values.is_null() {
            if irow.is_null() || jcol.is_null() {
                return cipopt::FALSE;
            }
            let mut k: usize = 0;
            for (row, col, _val) in p.hphi().iter() {
                *irow.add(k) = row.to_i32().unwrap();
                *jcol.add(k) = col.to_i32().unwrap();
                k += 1;
            }
            for (row, col, _val) in p.hcomb().iter() {
                *irow.add(k) = row.to_i32().unwrap();
                *jcol.add(k) = col.to_i32().unwrap();
                k += 1;
            }
        }
        else {
            if new_x == cipopt::TRUE {
                if x.is_null() {
                    return cipopt::FALSE;
                }
                p.evaluate(slice::from_raw_parts(x, p.nx()));
            }
            if new_lambda == cipopt::TRUE {
                if lambda.is_null() {
                    return cipopt::FALSE;
                }
                p.combine_h(slice::from_raw_parts(lambda.add(p.na()), p.nf()));
            }
            if new_x == cipopt::TRUE || new_lambda == cipopt::TRUE {
                let mut k: usize = 0;
                for (_row, _col, val) in p.hphi().iter() {
                    *values.add(k) = obj_factor*(*val); 
                    k += 1;
                }
                for (_row, _col, val) in p.hcomb().iter() {
                    *values.add(k) = *val; 
                    k += 1;
                }
            }   
        }
    };
    cipopt::TRUE
}

unsafe fn add_int_option(cprob: cipopt::IpoptProblem, 
                         key: &str, 
                         val: i32) 
                         -> Result<(), SimpleError> {
    let cstr = CString::new(key).unwrap();
    match cipopt::AddIpoptIntOption(cprob, 
                                    cstr.as_ptr() as *const c_char, 
                                    val) {
        cipopt::TRUE => Ok(()),
        _ =>  return Err(SimpleError::new(format!("unable to set parameter {}", key)))
    }
}

unsafe fn add_str_option(cprob: cipopt::IpoptProblem, 
                         key: &str, 
                         val: &str) 
                         -> Result<(), SimpleError> {
    let ckey = CString::new(key).unwrap();
    let cval = CString::new(val).unwrap();
    match cipopt::AddIpoptStrOption(cprob, 
                                    ckey.as_ptr() as *const c_char, 
                                    cval.as_ptr() as *const c_char) {
        cipopt::TRUE => Ok(()),
        _ =>  return Err(SimpleError::new(format!("unable to set parameter {}", key)))
    }
}

#[cfg(test)]
mod tests {

    use serial_test::serial;

    use crate::matrix::coo::CooMat;
    use crate::problem::base::Problem;
    use crate::problem::lp::ProblemLp;
    use crate::problem::nlp::ProblemNlp;
    use crate::solver::base::{Solver, SolverParam, SolverStatus};
    use crate::solver::ipopt::SolverIpopt;
    use crate::assert_vec_approx_eq;

    #[test]
    #[serial]
    fn ipopt_solve_nlp() {

        // Sample NLP problem 
        // min        x0*x3*(x0+x1+x2) + x2 = x0*x3*x0 + x0*x3*x1 + x0*x3*x2 + x2
        // subject to x0*x1*x2*x3 - x4 == 0
        //            x0*x0 + x1*x1 + x2*x2 + x3*x3 - 40 == 0
        //             1 <= x0 <= 5
        //             1 <= x1 <= 5
        //             1 <= x2 <= 5
        //             1 <= x3 <= 5
        //            25 <= x4 <= inf

        // x0
        let x0 = vec![1., 5., 5., 1., 0.];

        // hphi
        // 2*x3         x3  x3 (2*x0+x1+x2) 0
        // x3           0   0  x0           0
        // x3           0   0  x0           0
        // (2*x0+x1+x2) x0  x0 0            0
        // 0            0   0  0            0
        let hphi = CooMat::from_pattern(
            (5, 5),
            vec![0, 1, 2, 3, 3, 3],
            vec![0, 0, 0, 0, 1, 2],
        );

        let a = CooMat::from_nnz((0, 5), 0);
        let b = Vec::new();

        // j
        // x1*x2*x3 x0*x2*x3 x0*x1*x3 x0*x1*x2 -1 
        // 2*x0     2*x1     2*x2     2*x3      0
        let j = CooMat::from_pattern(
            (2, 5),
            vec![0, 0, 0, 0, 0, 1, 1, 1, 1],
            vec![0, 1, 2, 3, 4, 0, 1, 2, 3],
        );

        // h0
        // 0     x2*x3 x1*x3 x1*x2 0 
        // x2*x3 0     x0*x3 x0*x2 0
        // x1*x3 x0*x3 0     x0*x1 0
        // x1*x2 x0*x2 x0*x1 0     0
        // 0     0     0     0     0
        // h1
        // 2 0 0 0 0 
        // 0 2 0 0 0
        // 0 0 2 0 0
        // 0 0 0 2 0
        // 0 0 0 0 0
        let h = vec![
            CooMat::from_pattern(
                (5, 5),
                vec![1, 2, 2, 3, 3, 3],
                vec![0, 0, 1, 0, 1, 2],
            ),
            CooMat::from_pattern(
                (5, 5),
                vec![0, 1, 2, 3],
                vec![0, 1, 2, 3],
            )
        ];

        // l
        let l = vec![1., 1., 1., 1., 25.];

        // u
        let u = vec![5., 5., 5., 5., 1e8];
        
        // eval_fn
        let eval_fn = Box::new(move | phi: &mut f64, 
                                      gphi: &mut Vec<f64>, 
                                      hphi: &mut CooMat<f64>,
                                      f: &mut Vec<f64>,
                                      j: &mut CooMat<f64>,
                                      h: &mut Vec<CooMat<f64>>,
                                      x: &[f64] | {

            assert_eq!(gphi.len(), x.len());

            let x0 = x[0];
            let x1 = x[1];
            let x2 = x[2];
            let x3 = x[3];
            let x4 = x[4];

            // phi
            *phi = x0*x3*(x0+x1+x2) + x2;

            // gphi
            gphi[0] = 2.*x0*x3 + x1*x3 + x2*x3;
            gphi[1] = x0*x3;
            gphi[2] = x0*x3 + 1.;
            gphi[3] = x0*(x0+x1+x2);
            gphi[4] = 0.;

            // hphi
            let hphi_data = hphi.data_mut();
            hphi_data[0] = 2.*x3;       // x0, x0
            hphi_data[1] = x3;          // x1, x0
            hphi_data[2] = x3;          // x2, x0
            hphi_data[3] = 2.*x0+x1+x2; // x3, x0
            hphi_data[4] = x0;          // x3, x1
            hphi_data[5] = x0;          // x3, x2

            // f
            f[0] = x0*x1*x2*x3 - x4;
            f[1] = x0*x0 + x1*x1 + x2*x2 + x3*x3 - 40.;

            // j
            let j_data = j.data_mut();
            j_data[0] = x1*x2*x3; // 0, x0
            j_data[1] = x0*x2*x3; // 0, x1
            j_data[2] = x0*x1*x3; // 0, x2
            j_data[3] = x0*x1*x2; // 0, x3
            j_data[4] = -1.;      // 0, x4
            j_data[5] = 2.*x0;    // 1, x0
            j_data[6] = 2.*x1;    // 1, x1
            j_data[7] = 2.*x2;    // 1, x2
            j_data[8] = 2.*x3;    // 1, x3

            // h0
            let h0_data = h[0].data_mut();
            h0_data[0] = x2*x3;
            h0_data[1] = x1*x3;
            h0_data[2] = x0*x3;
            h0_data[3] = x1*x2;
            h0_data[4] = x0*x2;
            h0_data[5] = x0*x1;

            // h1
            let h1_data = h[1].data_mut();
            h1_data[0] = 2.;
            h1_data[1] = 2.;
            h1_data[2] = 2.;
            h1_data[3] = 2.;
        });

        let mut p = Problem::Nlp(ProblemNlp::new(
            hphi, 
            a,
            b,
            j,
            h,
            l,
            u,
            Some(x0),
            eval_fn
        ));

        let mut s = SolverIpopt::new();
        s.set_param("print_level", SolverParam::IntParam(0)).unwrap();
        s.set_param("sb", SolverParam::StrParam("yes".to_string())).unwrap();
        let (status, solution) = s.solve(&mut p).unwrap();

        assert_eq!(status, SolverStatus::Solved);
        assert_vec_approx_eq!(solution.x,
                              &vec![1., 4.742999629, 3.821149993, 1.379408294, 25.],
                              epsilon=1e-7);

    }

    #[test]
    #[serial]
    fn ipopt_solve_lp() {

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
 
        let mut s = SolverIpopt::new();
        s.set_param("print_level", SolverParam::IntParam(0)).unwrap();
        s.set_param("sb", SolverParam::StrParam("yes".to_string())).unwrap();
        let (status, solution) = s.solve(&mut p).unwrap();

        assert_eq!(status, SolverStatus::Solved);
        assert_vec_approx_eq!(solution.x, 
                              &vec![1.7142857, 2.8571429, -1.1428571, 0., 0.], 
                              epsilon=1e-6);
        assert_vec_approx_eq!(solution.lam, 
                              &vec![0., 31.428571, 21.428571], 
                              epsilon=1e-6);
        assert_vec_approx_eq!(solution.mu, 
                              &vec![1.4210855e-14, 0., 0., 3.1428571e+01, 2.1428571e+01], 
                              epsilon=1e-6);
        assert_vec_approx_eq!(solution.pi, 
                              &vec![0.;5], 
                              epsilon=1e-6);
    }
}