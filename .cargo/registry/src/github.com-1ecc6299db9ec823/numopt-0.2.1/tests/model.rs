use serial_test::serial;

use numopt::model::*;
use numopt::solver::*;

#[test]
#[serial]
fn numopt_model_solve_milp_cbc_cmd() {

    let x1 = VariableScalar::new_integer("x1");
    let x2 = VariableScalar::new_integer("x2");
    let x3 = VariableScalar::new_continuous("x3");
    let x4 = VariableScalar::new_continuous("x4");

    let f = -&x1 - &x2;
    let c1 = (-2.*&x1+ 2.*&x2 + &x3).equal(1.);
    let c2 = (-8.*&x1 + 10.*&x2 + &x4).equal(13.);
    let c3 = &x4.geq(0.);
    let c4 = &x3.leq(0.);

    let mut m = Model::new();
    m.set_objective(Objective::minimize(&f));
    m.add_constraint(&c1);
    m.add_constraint(&c2);
    m.add_constraint(&c3);
    m.add_constraint(&c4);

    let mut s = SolverCbcCmd::new();
    s.set_param("logLevel", SolverParam::IntParam(0)).unwrap();
    m.solve(&s).unwrap();

    assert_eq!(*m.solver_status().unwrap(), SolverStatus::Solved);
    assert_eq!(*m.final_primals().get(&x1).unwrap(), 1.);
    assert_eq!(*m.final_primals().get(&x2).unwrap(), 2.);
}

#[cfg(feature = "ipopt")] 
#[test]
#[serial]
fn numopt_model_solve_nlp1_ipopt() {

    // Rosenbrock

    use approx::assert_abs_diff_eq;

    const N: usize = 100;

    let x: Vec<Node>= (0..N).map(|i| VariableScalar::new_continuous(format!("x{}", i).as_ref()))
                            .collect();

    let mut f = ConstantScalar::zero();
    for i in 0..(N-1) {
        f = f + 100.*(&x[i+1] - &x[i]*&x[i])*(&x[i+1] - &x[i]*&x[i]) + (1. - &x[i])*(1. - &x[i]);
    }

    let mut m = Model::new();
    m.set_objective(Objective::minimize(&f));

    let mut s = SolverIpopt::new();
    s.set_param("print_level", SolverParam::IntParam(0)).unwrap();
    s.set_param("sb", SolverParam::StrParam("yes".to_string())).unwrap();
    m.solve(&s).unwrap();

    assert_eq!(*m.solver_status().unwrap(), SolverStatus::Solved);
    let final_primals = m.final_primals();
    for xx in x.iter() {
        assert_abs_diff_eq!(*final_primals.get(xx).unwrap(), 1., epsilon=1e-8);
    }
}

#[cfg(feature = "ipopt")] 
#[test]
#[serial]
fn numopt_model_solve_nlp2_ipopt() {

    // https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/clnlbeam.jl

    use num_traits::ToPrimitive;
    use approx::assert_abs_diff_eq;

    const N: usize = 100;
    const ALPHA: f64 = 350.;
    let h: f64 = 1./N.to_f64().unwrap();

    let t: Vec<Node>= (0..N+1).map(|i| VariableScalar::new_continuous(format!("t{}", i).as_ref()))
                              .collect();    
    let x: Vec<Node>= (0..N+1).map(|i| VariableScalar::new_continuous(format!("x{}", i).as_ref()))
                              .collect();
    let u: Vec<Node>= (0..N+1).map(|i| VariableScalar::new_continuous(format!("u{}", i).as_ref()))
                              .collect();

    let mut f = ConstantScalar::zero();
    for i in 0..N {
        f = f + 0.5*h*(&u[i]*&u[i] + &u[i+1]*&u[i+1]) +
                0.5*ALPHA*h*(t[i].cos() + t[i+1].cos());
    }
    
    let mut m = Model::new();
    m.set_objective(Objective::minimize(&f));
    for i in 0..N {
        m.add_constraint(&(&x[i+1] - &x[i] - 0.5*h*(t[i+1].sin() + t[i].sin())).equal(0.));
        m.add_constraint(&(&t[i+1] - &t[i] - 0.5*h*(&u[i+1] - &u[i])).equal(0.));
    }
    for i in 0..(N+1) {
        m.add_constraint(&t[i].leq(1.));
        m.add_constraint(&t[i].geq(-1.));
        m.add_constraint(&x[i].leq(0.05));
        m.add_constraint(&x[i].geq(-0.05));
    }

    let mut s = SolverIpopt::new();
    s.set_param("print_level", SolverParam::IntParam(0)).unwrap();
    s.set_param("sb", SolverParam::StrParam("yes".to_string())).unwrap();
    m.solve(&s).unwrap();

    assert_eq!(*m.solver_status().unwrap(), SolverStatus::Solved);
    assert_abs_diff_eq!(f.evaluate(&m.final_primals()), 350., epsilon=1e-8);
}
