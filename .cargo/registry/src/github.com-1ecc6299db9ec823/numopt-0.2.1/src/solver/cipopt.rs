#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use libc::{c_int, 
           c_void, 
           c_char,
           c_double};

#[repr(C)] pub struct IpoptProblemInfo { _private: [u8; 0] }

pub const TRUE: c_int = 1;
pub const FALSE: c_int = 0;

pub type IpoptProblem = *mut IpoptProblemInfo;

type Eval_F_CB = extern fn(c_int, 
                           *const c_double, 
                           c_int, 
                           *mut c_double, 
                           *mut c_void) -> c_int;

type Eval_Grad_F_CB = extern fn(c_int, 
                                *const c_double, 
                                c_int,
                                *mut c_double, 
                                *mut c_void) -> c_int;

type Eval_G_CB = extern fn(c_int, 
                           *const c_double, 
                           c_int,
                           c_int, 
                           *mut c_double, 
                           *mut c_void) -> c_int;

type Eval_Jac_G_CB = extern fn(c_int,
                               *const c_double, 
                               c_int,
                               c_int, 
                               c_int,
                               *mut c_int, 
                               *mut c_int, 
                               *mut c_double,
                               *mut c_void) -> c_int;

type Eval_H_CB = extern fn(c_int, 
                           *const c_double, 
                           c_int, 
                           c_double,
                           c_int, 
                           *const c_double, 
                           c_int,
                           c_int,
                           *mut c_int, 
                           *mut c_int, 
                           *mut c_double,
                           *mut c_void) -> c_int;

#[link(name = "ipopt")]
extern {

    pub fn CreateIpoptProblem(n: c_int,
                              x_L: *const c_double,
                              x_U: *const c_double,
                              m: c_int,
                              g_L: *const c_double,
                              g_U: *const c_double,
                              nele_jac: c_int,
                              nele_hess: c_int,
                              index_style: c_int,
                              eval_f: Eval_F_CB,
                              eval_g: Eval_G_CB,
                              eval_grad_f: Eval_Grad_F_CB,
                              eval_jac_g: Eval_Jac_G_CB,
                              eval_h: Eval_H_CB
                             ) -> IpoptProblem;

    pub fn FreeIpoptProblem(ipopt_problem: IpoptProblem) -> ();

    pub fn IpoptSolve(ipopt_problem: IpoptProblem,
                      x: *mut c_double,
                      g: *mut c_double,
                      obj_value: *mut c_double,
                      mult_g: *mut c_double,
                      mult_x_L: *mut c_double,
                      mult_x_U: *mut c_double,
                      user_data: *mut c_void) -> c_int;

    pub fn AddIpoptIntOption(ipopt_problem: IpoptProblem, 
                             keyword: *const c_char, 
                             val: c_int) -> c_int;

    pub fn AddIpoptStrOption(ipopt_problem: IpoptProblem, 
                             keyword: *const c_char, 
                             val: *const c_char) -> c_int;
}


