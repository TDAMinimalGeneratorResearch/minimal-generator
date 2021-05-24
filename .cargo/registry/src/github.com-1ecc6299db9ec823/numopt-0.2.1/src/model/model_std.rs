//! Structures and traits for transforming optimization models to standard form.

use std::collections::{HashSet, HashMap};

use crate::matrix::coo::CooMat;
use crate::problem::base::Problem;
use crate::problem::minlp::ProblemMinlp;
use crate::problem::milp::ProblemMilp;
use crate::problem::nlp::ProblemNlp;
use crate::problem::lp::ProblemLp;

use crate::model::node::Node;
use crate::model::node_base::NodeBase;
use crate::model::node_std::{NodeStd, NodeStdComp};
use crate::model::constant::ConstantScalar;
use crate::model::constraint::Constraint;
use crate::model::constraint_std::{ConstraintStd, ConstraintStdComp};
use crate::model::model::{Model, Objective};

const INF: f64 = 1e8;

/// Optimization model standard components.
pub struct ModelStdComp {

    /// Standard components of objective expression.
    pub obj: NodeStdComp,

    /// Standard components of constraints.
    pub constr: ConstraintStdComp,
}

/// Optimization model standard problem.
pub struct ModelStdProb {

    /// Problem in standard form.
    pub prob: Problem,

    /// Map between optimization model variables and their index
    /// in the vector of variables for the problem in standard form.
    pub var2index: HashMap<Node, usize>,

    /// Map between linear equality constraint row and model constraint.
    pub aindex2constr: HashMap<usize, Constraint>,
    
    /// Map between nonlinear equality constraint row and model constraint.
    pub jindex2constr: HashMap<usize, Constraint>,

    /// Map between variable upper limit and model constraint.
    pub uindex2constr: HashMap<usize, Constraint>,

    /// Map between variable lower limit and model constraint.
    pub lindex2constr: HashMap<usize, Constraint>,
}

/// A trait for transforming optimization models to problems in standard form.
pub trait ModelStd {

    /// Obtains standard components of model.
    fn std_components(&self) -> ModelStdComp;

    /// Obtains standard problem of model.
    fn std_problem(&self) -> ModelStdProb;
}

impl ModelStd for Model {

    fn std_components(&self) -> ModelStdComp {

        // Objective std comp
        let obj = match self.objective() {
            Objective::Maximize(f) => (-f).std_components(),
            Objective::Minimize(f) => f.std_components(),
            Objective::Empty => ConstantScalar::new(0.).std_components(),
        };

        // Constraint std comp
        let mut arow: usize = 0;
        let mut jrow: usize = 0;
        let mut constr = ConstraintStdComp::new();
        for c in self.constraints().iter() {
            constr += c.std_components(&mut arow, &mut jrow);
        }

        // Return
        ModelStdComp {
            obj: obj,
            constr: constr,
        }
    }

    fn std_problem(&self) -> ModelStdProb {

        // Components
        let comp = self.std_components();

        // Variables
        let mut varset: HashSet<Node> = comp.obj.prop.a.keys()
                                                       .map(|x| x.clone())
                                                       .collect();
        for p in comp.constr.prop.iter() {
            varset.extend(p.a.keys().map(|x| x.clone()));
        }
        let num_vars: usize = varset.len();
        let mut vars: Vec<Node> = varset.into_iter().collect();
        vars.sort_by(|x,y| x.name().partial_cmp(y.name()).unwrap());
        let var2index: HashMap<Node, usize> = vars.into_iter()
                                                  .enumerate()
                                                  .map(|(i,v)| (v,i))
                                                  .collect();
        let var2index_eval: HashMap<Node, usize> = var2index.iter()
                                                            .map(|(v,i)| (v.clone(), *i))
                                                            .collect();
    
        // Objective (phi)
        let phi_data = comp.obj.phi;
        let mut gphi_indices: Vec<usize> = Vec::with_capacity(comp.obj.gphi.len());
        let mut gphi_data: Vec<Node> = Vec::with_capacity(comp.obj.gphi.len());
        for (v, e) in comp.obj.gphi.into_iter() {
            gphi_indices.push(*var2index.get(&v).unwrap());
            gphi_data.push(e);
        }
        let mut hphi_row: Vec<usize> = Vec::with_capacity(comp.obj.hphi.len());
        let mut hphi_col: Vec<usize> = Vec::with_capacity(comp.obj.hphi.len());
        let mut hphi_data: Vec<Node> = Vec::with_capacity(comp.obj.hphi.len());
        for (v1, v2, e) in comp.obj.hphi.into_iter() {
            let i = *var2index.get(&v1).unwrap();
            let j = *var2index.get(&v2).unwrap();
            if i >= j {
                hphi_row.push(i);
                hphi_col.push(j);
            }
            else {
                hphi_row.push(j);
                hphi_col.push(i);
            }
            hphi_data.push(e);
        }
        let hphi_mat = CooMat::new(
            (num_vars, num_vars),
            hphi_row,
            hphi_col,
            vec![0.; hphi_data.len()]
        );

        // Objective grad (c)
        let mut c_data: Vec<f64> = vec![0.; num_vars];
        for (var, val) in comp.obj.prop.a.iter() {
            c_data[*var2index.get(var).unwrap()] = *val;
        }

        // Linear equality constraints (Ax = b)
        let aindex2constr: HashMap<usize, Constraint> = comp.constr.ca.into_iter()
                                                                      .enumerate()
                                                                      .collect();
        let num_a: usize = comp.constr.b.len();
        let mut a_row: Vec<usize> = Vec::with_capacity(comp.constr.a.len());
        let mut a_col: Vec<usize> = Vec::with_capacity(comp.constr.a.len());
        let mut a_data: Vec<f64> = Vec::with_capacity(comp.constr.a.len());
        for (row, var, val) in comp.constr.a.into_iter() {
            a_row.push(row);
            a_col.push(*var2index.get(&var).unwrap());
            a_data.push(val);
        }
        let a_mat = CooMat::new(
            (num_a, num_vars),
            a_row,
            a_col,
            a_data
        );

        let b_data = comp.constr.b;

        // Nonlinear equality constraints (f(x) = 0)
        let jindex2constr: HashMap<usize, Constraint> = comp.constr.cj.into_iter()
                                                                      .enumerate()
                                                                      .collect();
        let num_j: usize = comp.constr.f.len();
        let mut j_row: Vec<usize> = Vec::with_capacity(comp.constr.j.len());
        let mut j_col: Vec<usize> = Vec::with_capacity(comp.constr.j.len());
        let mut j_data: Vec<Node> = Vec::with_capacity(comp.constr.j.len());
        for (row, var, exp) in comp.constr.j.into_iter() {
            j_row.push(row);
            j_col.push(*var2index.get(&var).unwrap());
            j_data.push(exp);
        }
        let j_mat = CooMat::new(
            (num_j, num_vars),
            j_row,
            j_col,
            vec![0.; j_data.len()]
        );
        let f_data = comp.constr.f;
        let mut h_data: Vec<Vec<Node>> = Vec::with_capacity(num_j);
        let mut h_vec: Vec<CooMat<f64>> = Vec::with_capacity(num_j);
        for hh in comp.constr.h.into_iter() {
            let mut hh_row: Vec<usize> = Vec::with_capacity(hh.len());
            let mut hh_col: Vec<usize> = Vec::with_capacity(hh.len());
            let mut hh_data: Vec<Node> = Vec::with_capacity(hh.len());
            for (v1, v2, exp) in hh.into_iter() {
                let i = *var2index.get(&v1).unwrap();
                let j = *var2index.get(&v2).unwrap();
                if i >= j {
                    hh_row.push(i);
                    hh_col.push(j);
                }
                else {
                    hh_row.push(j);
                    hh_col.push(i);
                }
                hh_data.push(exp);
            }
            h_vec.push(CooMat::new(
                (num_vars, num_vars),
                hh_row,
                hh_col,
                vec![0.; hh_data.len()]
            ));
            h_data.push(hh_data);
        }

        // Bounds (l <= x <= u)
        let mut uindex2constr: HashMap<usize, Constraint> = HashMap::new();
        let mut lindex2constr: HashMap<usize, Constraint> = HashMap::new();
        let mut u_data = vec![INF; num_vars];
        let mut l_data = vec![-INF; num_vars];
        for (var, val, constr) in comp.constr.u.into_iter() {
            let index = *var2index.get(&var).unwrap();
            if val <= u_data[index] {
                u_data[index] = val;
                uindex2constr.insert(index, constr);
            }
        }
        for (var, val, constr) in comp.constr.l.into_iter() {
            let index = *var2index.get(&var).unwrap();
            if val >= l_data[index] {
                l_data[index] = val;
                lindex2constr.insert(index, constr);
            }
        }

        // Integer restrictions
        let mut num_int: usize = 0;
        let mut p_data = vec![false; num_vars];
        for (var, index) in var2index.iter() {
            match var {
                Node::VariableScalar(x) => {
                    if x.is_integer() {
                        p_data[*index] = true;
                        num_int += 1;
                    }
                }
                _ => (),
            }
        }

        // Initial values
        let mut x0_data: Vec<f64> = vec![0.; num_vars];
        for (var, val) in self.init_primals().iter() {
            match var2index.get(var) {
                Some(index) => x0_data[*index] = *val,
                None => (), 
            }
        }
       
        // Eval
        let eval_fn = Box::new(move | phi: &mut f64, 
                                      gphi: &mut Vec<f64>, 
                                      hphi: &mut CooMat<f64>,
                                      f: &mut Vec<f64>,
                                      j: &mut CooMat<f64>,
                                      h: &mut Vec<CooMat<f64>>,
                                      x: &[f64] | {

            // Var values
            let mut var_values: HashMap<&Node, f64> = HashMap::with_capacity(x.len());
            for (var, index) in var2index_eval.iter() {
                var_values.insert(var, x[*index]);
            }

            // phi
            *phi = phi_data.evaluate(&var_values);

            // gphi
            for (index, exp) in gphi_indices.iter().zip(gphi_data.iter()) {
                (*gphi)[*index] = exp.evaluate(&var_values);
            }

            // hphi
            let hphi_dest = hphi.data_mut();
            for (val, exp) in hphi_dest.iter_mut().zip(hphi_data.iter()) {
                *val = exp.evaluate(&var_values);           
            }

            // f
            for (val, exp) in f.iter_mut().zip(f_data.iter()) {
                *val = exp.evaluate(&var_values);
            }

            // j
            let j_dest = j.data_mut();
            for (val, exp) in j_dest.iter_mut().zip(j_data.iter()) {
                *val = exp.evaluate(&var_values);
            }
            
            // h
            for (hh, hh_data) in h.iter_mut().zip(h_data.iter()) {
                let hh_dest = hh.data_mut();
                for (val, exp) in hh_dest.iter_mut().zip(hh_data.iter()) {
                    *val = exp.evaluate(&var_values);
                }
            }
        });

        // Problem
        let problem: Problem;

        // Lp
        if comp.obj.prop.affine && num_j == 0 && num_int == 0 {
            problem = Problem::Lp(
                ProblemLp::new(
                    c_data,
                    a_mat,
                    b_data,
                    l_data,
                    u_data,
                    Some(x0_data)
                )
            );
        }

        // Milp
        else if comp.obj.prop.affine && num_j == 0 && num_int > 0 {
            problem = Problem::Milp(
                ProblemMilp::new(
                    c_data,
                    a_mat,
                    b_data,
                    l_data,
                    u_data,
                    p_data,
                    Some(x0_data)
                )
            );
        }

        // Nlp
        else if num_int == 0 {
            problem = Problem::Nlp(
                ProblemNlp::new(
                    hphi_mat,
                    a_mat,
                    b_data,
                    j_mat,
                    h_vec,
                    l_data,
                    u_data,
                    Some(x0_data),
                    eval_fn,
                )
            );
        }

        // Minlp
        else {
            problem = Problem::Minlp(
                ProblemMinlp::new(
                    hphi_mat,
                    a_mat,
                    b_data,
                    j_mat,
                    h_vec,
                    l_data,
                    u_data,
                    p_data,
                    Some(x0_data),
                    eval_fn,
                )
            );
        }

        // Return
        ModelStdProb {
            prob: problem,
            var2index: var2index,
            aindex2constr: aindex2constr,
            jindex2constr: jindex2constr,
            uindex2constr: uindex2constr,
            lindex2constr: lindex2constr,
        }
    }
}

#[cfg(test)]
mod tests {

    use maplit::hashmap;

    use super::*;
    use crate::problem::base::Problem;
    use crate::model::node_func::NodeFunc;
    use crate::model::node_cmp::NodeCmp;
    use crate::model::variable::VariableScalar;
    use crate::assert_vec_approx_eq;

    #[test]
    fn model_std_problem_lp() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let c1 = (2.*&x + &y).equal(2.);
        let c2 = x.leq(5.);
        let c3 = x.geq(0.);
        let c4 = y.leq(5.);
        let c5 = y.geq(0.);
        let c6 = (7.*&x + 5.*&y - 10.).geq(&x + &y - 3.);

        let mut m = Model::new();
        m.set_objective(Objective::maximize(&(3.*&x + 4.*&y + 1.)));
        m.add_constraint(&c1);
        m.add_constraint(&c2);
        m.add_constraint(&c3);
        m.add_constraint(&c4);
        m.add_constraint(&c5);
        m.add_constraint(&c6);
        m.set_init_primals(&hashmap!{ &x => 2., &y => 3. });

        let std_p = m.std_problem();
        let lp = match std_p.prob {
            Problem::Lp(x) => x,
            _ => panic!("invalid std problem")
        };

        assert_vec_approx_eq!(lp.x0().unwrap(), vec![0., 2., 3.], epsilon=0.);
        assert_vec_approx_eq!(lp.c(), vec![0., -3., -4.,], epsilon=0.);
        assert_eq!(lp.na(), 2);
        assert_eq!(lp.nx(), 3);
        assert_eq!(lp.a().nnz(), 5);
        for (row, col, val) in lp.a().iter() {
            if *row == 0 && *col == 1 {
                assert_eq!(*val, 2.);
            }
            else if *row == 0 && *col == 2 {
                assert_eq!(*val, 1.);
            }
            else if *row == 1 && *col == 0 {
                assert_eq!(*val, -1.);
            }
            else if *row == 1 && *col == 1 {
                assert_eq!(*val, 6.);
            }
            else if *row == 1 && *col == 2 {
                assert_eq!(*val, 4.);
            }
            else {
                panic!("invalid a matrix")
            }
        }
        assert_vec_approx_eq!(lp.b(), vec![2., 7.], epsilon=0.);
        assert_vec_approx_eq!(lp.l(), vec![0., 0., 0.], epsilon=0.);
        assert_vec_approx_eq!(lp.u(), vec![1e8, 5., 5.], epsilon=0.);

        assert_eq!(std_p.var2index.len(), 3);
        for (var, index) in std_p.var2index.iter() {
            if *var == x {
                assert_eq!(*index, 1);
            }
            else if *var == y {
                assert_eq!(*index, 2);
            }
            else if (*var).name() == "_s_a1_" {
                assert_eq!(*index, 0);
            }
            else {
                panic!("invalid var2index data");
            }
        }
        assert_eq!(std_p.aindex2constr.len(), 2);
        assert_eq!(*std_p.aindex2constr.get(&0).unwrap(), c1);
        assert_eq!(*std_p.aindex2constr.get(&1).unwrap(), c6);
        assert_eq!(std_p.jindex2constr.len(), 0);
        assert_eq!(std_p.uindex2constr.len(), 2);
        assert_eq!(*std_p.uindex2constr.get(&1).unwrap(), c2);
        assert_eq!(*std_p.uindex2constr.get(&2).unwrap(), c4);
        assert_eq!(std_p.lindex2constr.len(), 3);
        assert_eq!(*std_p.lindex2constr.get(&0).unwrap(), c6);
        assert_eq!(*std_p.lindex2constr.get(&1).unwrap(), c3);
        assert_eq!(*std_p.lindex2constr.get(&2).unwrap(), c5);
    }

    #[test]
    fn model_std_problem_milp() {

        let x = VariableScalar::new_integer("x");
        let y = VariableScalar::new_continuous("y");

        let c1 = (2.*&x + &y).equal(2.);
        let c2 = x.leq(5.);
        let c3 = x.geq(0.);
        let c4 = y.leq(5.);
        let c5 = y.geq(0.);

        let mut m = Model::new();
        m.set_objective(Objective::minimize(&(3.*&x + 4.*&y + 1.)));
        m.add_constraint(&c1);
        m.add_constraint(&c2);
        m.add_constraint(&c3);
        m.add_constraint(&c4);
        m.add_constraint(&c5);
        m.set_init_primals(&hashmap!{ &x => 2., &y => 3. });

        let std_p = m.std_problem();
        let milp = match std_p.prob {
            Problem::Milp(x) => x,
            _ => panic!("invalid std problem")
        };

        assert_vec_approx_eq!(milp.x0().unwrap(), vec![2., 3.], epsilon=0.);
        assert_vec_approx_eq!(milp.c(), vec![3., 4.,], epsilon=0.);
        assert_eq!(milp.na(), 1);
        assert_eq!(milp.nx(), 2);
        assert_eq!(milp.a().nnz(), 2);
        for (row, col, val) in milp.a().iter() {
            if *row == 0 && *col == 0 {
                assert_eq!(*val, 2.);
            }
            else if *row == 0 && *col == 1 {
                assert_eq!(*val, 1.);
            }
            else {
                panic!("invalid a matrix")
            }
        }
        assert_vec_approx_eq!(milp.b(), vec![2.], epsilon=0.);
        assert_vec_approx_eq!(milp.l(), vec![0., 0.], epsilon=0.);
        assert_vec_approx_eq!(milp.u(), vec![5., 5.], epsilon=0.);
        assert_eq!(milp.p().len(), 2);
        assert!(milp.p()[0]);
        assert!(!milp.p()[1]);

        assert_eq!(std_p.var2index.len(), 2);
        assert_eq!(*std_p.var2index.get(&x).unwrap(), 0);
        assert_eq!(*std_p.var2index.get(&y).unwrap(), 1);
        assert_eq!(std_p.aindex2constr.len(), 1);
        assert_eq!(*std_p.aindex2constr.get(&0).unwrap(), c1);
        assert_eq!(std_p.jindex2constr.len(), 0);
        assert_eq!(std_p.uindex2constr.len(), 2);
        assert_eq!(*std_p.uindex2constr.get(&0).unwrap(), c2);
        assert_eq!(*std_p.uindex2constr.get(&1).unwrap(), c4);
        assert_eq!(std_p.lindex2constr.len(), 2);
        assert_eq!(*std_p.lindex2constr.get(&0).unwrap(), c3);
        assert_eq!(*std_p.lindex2constr.get(&1).unwrap(), c5);
    }

    #[test]
    fn model_std_problem_nlp() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let f = 3.*x.cos() + 4.*&y*&y - 10.;
        let c1 = (2.*&x + 3.*&y).equal(7.);
        let c2 = (x.sin() + &x*&y + 4.).leq(&x + 10.);
        let c3 = 5_f64.geq(&x*&x);

        let mut m = Model::new();
        m.set_objective(Objective::maximize(&f));
        m.add_constraint(&c1);
        m.add_constraint(&c2);
        m.add_constraint(&c3);
        m.set_init_primals(&hashmap!{ &x => 2., &y => 3. });

        let std_p = m.std_problem();
        let mut nlp = match std_p.prob {
            Problem::Nlp(x) => x,
            _ => panic!("invalid std problem")
        };

        nlp.evaluate(&vec![2., 3., 4., 5.]);
        nlp.combine_h(&vec![1.5, 2.3]);

        assert_vec_approx_eq!(nlp.x0().unwrap(), vec![0., 0., 2., 3.], epsilon=0.);
        assert_eq!(nlp.phi(), -(3.*4_f64.cos() + 4.*5.*5. - 10.));
        assert_vec_approx_eq!(nlp.gphi(), 
                              vec![0., 0., 3.*4_f64.sin(), -8.*5.],
                              epsilon=1e-8);
        assert_eq!(nlp.hphi().nnz(), 2);
        for (v1, v2, val) in nlp.hphi().iter() {
            if *v1 == 2 && *v2 == 2 {
                assert_eq!(*val, 3.*4_f64.cos());
            }
            else if *v1 == 3 && *v2 == 3 {
                assert_eq!(*val, -8.);
            }
            else {
                panic!("invalid variable pair");
            }
        }

        assert_eq!(nlp.a().nnz(), 2);
        for (row, col, val) in nlp.a().iter() {
            if *row == 0 && *col == 2 {
                assert_eq!(*val, 2.);
            }
            else if *row == 0 && *col == 3 {
                assert_eq!(*val, 3.);
            }
            else {
                panic!("invalid a matrix entry");
            }
        }
        assert_vec_approx_eq!(nlp.b(), vec![7.], epsilon=0.);

        assert_eq!(nlp.f().len(), 2);
        assert_eq!(nlp.f()[0], 4_f64.sin() + 4.*5. + 4. - 4. - 10. - 2.);
        assert_eq!(nlp.f()[1], 5. - 4.*4. - 3.);
        assert_eq!(nlp.j().nnz(), 5);
        for (row, col, val) in nlp.j().iter() {
            if *row == 0 && *col == 2 {
                assert_eq!(*val, 4_f64.cos() + 5. - 1.);
            }
            else if *row == 0 && *col == 3 {
                assert_eq!(*val, 4.);
            }
            else if *row == 0 && *col == 0 {
                assert_eq!(*val, -1.);
            }
            else if *row == 1 && *col == 2 {
                assert_eq!(*val, -2.*4.);
            }
            else if *row == 1 && *col == 1 {
                assert_eq!(*val, -1.);
            }
            else {
                panic!("invalid j matrix entry");
            }
        }
        assert_eq!(nlp.h().len(), 2);
        assert_eq!(nlp.h()[0].nnz(), 2);
        assert_eq!(nlp.h()[1].nnz(), 1);
        for (var1, var2, val) in nlp.h()[0].iter() {
            if *var1 == 2 && *var2 == 2 {
                assert_eq!(*val, -4_f64.sin());
            }
            else if *var1 == 3 && *var2 == 2 {
                assert_eq!(*val, 1.);
            }
            else {
                panic!("invalid h[0] entry")
            }
        }
        for (var1, var2, val) in nlp.h()[1].iter() {
            if *var1 == 2 && *var2 == 2 {
                assert_eq!(*val, -2.);
            }
            else {
                panic!("invalid h[1] entry")
            }
        }
        assert_eq!(nlp.hcomb().nnz(), 3);
        let mut first = true;
        for (var1, var2, val) in nlp.hcomb().iter() {
            if *var1 == 2 && *var2 == 2 && first {
                assert_eq!(*val, 1.5*(-4_f64.sin()));
                first = false;
            }
            else if *var1 == 2 && *var2 == 2 && !first {
                assert_eq!(*val, 2.3*-2.);
            }
            else if *var1 == 3 && *var2 == 2 {
                assert_eq!(*val, 1.*1.5);
            }
            else {
                panic!("invalid hcomb entry")
            }
        }
        assert_vec_approx_eq!(nlp.l(), vec![-1e8, 0., -1e8, -1e8], epsilon=0.);
        assert_vec_approx_eq!(nlp.u(), vec![0., 1e8, 1e8, 1e8], epsilon=0.);

        assert_eq!(std_p.var2index.len(), 4);
        for (var, index) in std_p.var2index.iter() {
            if *var == x{
                assert_eq!(*index, 2);
            }
            else if *var == y {
                assert_eq!(*index, 3);
            }
            else if (*var).name() == "_s_j0_" {
                assert_eq!(*index, 0);
            }
            else if (*var).name() == "_s_j1_" {
                assert_eq!(*index, 1);
            }
            else {
                panic!("invalid var2index entry");
            }
        }
        assert_eq!(std_p.aindex2constr.len(), 1);
        assert_eq!(*std_p.aindex2constr.get(&0).unwrap(), c1);
        assert_eq!(std_p.jindex2constr.len(), 2);
        assert_eq!(*std_p.jindex2constr.get(&0).unwrap(), c2);
        assert_eq!(*std_p.jindex2constr.get(&1).unwrap(), c3);
        assert_eq!(std_p.uindex2constr.len(), 1);
        assert_eq!(*std_p.uindex2constr.get(&0).unwrap(), c2);
        assert_eq!(std_p.lindex2constr.len(), 1);
        assert_eq!(*std_p.lindex2constr.get(&1).unwrap(), c3);
    }
}