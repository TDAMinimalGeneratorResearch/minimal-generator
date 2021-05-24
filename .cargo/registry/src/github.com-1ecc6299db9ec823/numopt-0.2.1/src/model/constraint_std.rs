//! Structures and traits for transforming constraints to standard form.

use std::ops::AddAssign;

use crate::model::node::Node;
use crate::model::node_std::NodeStd;
use crate::model::node_std::NodeStdProp;
use crate::model::constant::ConstantScalar;
use crate::model::variable::VariableScalar;
use crate::model::constraint::{Constraint, ConstraintKind};

/// Constraint standard components.
pub struct ConstraintStdComp {

    /// List of constraints that map to linear equality constraints.
    pub ca: Vec<Constraint>,             // constraints

    /// List of constraints that map to nonlinear equality constraints.
    pub cj: Vec<Constraint>,             // constraints

    /// Data for Jacobian matrix for linear equality constraints.
    pub a: Vec<(usize, Node, f64)>,      // row, var, value

    /// Right-hand-side vector of linear equality constraint functions.
    pub b: Vec<f64>,                     // values

    /// Expressions for nonlinear equality constraints.
    pub f: Vec<Node>,                    // expressions

    /// Data for Jacobian of nonlinear equality constraint expressions.
    pub j: Vec<(usize, Node, Node)>,     // row, var, expression

    /// Data for Hessians of nonlinear equality constraint expressions.
    pub h: Vec<Vec<(Node, Node, Node)>>, // var, var, expression

    /// Data for variable upper limits.
    pub u: Vec<(Node, f64, Constraint)>, // var, value, constraint

    /// Data for variable lower limits.
    pub l: Vec<(Node, f64, Constraint)>, // var, value, constraint

    /// Standard properties of constraint expressions.
    pub prop: Vec<NodeStdProp>,
}

/// Trait for transforming constraints to standard form.
pub trait ConstraintStd {

    /// Gets standard components of constraint.
    fn std_components(&self, arow: &mut usize, jrow: &mut usize) -> ConstraintStdComp;
}

impl ConstraintStdComp {

    /// Creates new container of constraint standard components.
    pub fn new() -> ConstraintStdComp {
        ConstraintStdComp {
            ca: Vec::new(),
            cj: Vec::new(),
            a: Vec::new(),
            b: Vec::new(),
            f: Vec::new(),
            j: Vec::new(),
            h: Vec::new(),
            u: Vec::new(),
            l: Vec::new(),
            prop: Vec::new(),
        }
    }
}

impl AddAssign for ConstraintStdComp {
    fn add_assign(&mut self, other: Self) {
        self.ca.extend(other.ca);
        self.cj.extend(other.cj);
        self.a.extend(other.a);
        self.b.extend(other.b);
        self.f.extend(other.f);
        self.j.extend(other.j);
        self.h.extend(other.h);
        self.u.extend(other.u);
        self.l.extend(other.l);
        self.prop.extend(other.prop);
    }
}

impl ConstraintStd for Constraint {

    fn std_components(&self, arow: &mut usize, jrow: &mut usize) -> ConstraintStdComp {

        let mut ca: Vec<Constraint> = Vec::new();
        let mut cj: Vec<Constraint> = Vec::new();
        let mut a: Vec<(usize, Node, f64)> = Vec::new();
        let mut b: Vec<f64> = Vec::new();
        let mut f: Vec<Node> = Vec::new();    
        let mut j: Vec<(usize, Node, Node)> = Vec::new();
        let mut h: Vec<Vec<(Node, Node, Node)>> = Vec::new();
        let mut u: Vec<(Node, f64, Constraint)> = Vec::new();
        let mut l: Vec<(Node, f64, Constraint)> = Vec::new();
    
        let exp = self.lhs()-self.rhs();
        let comp = exp.std_components();
        let mut prop = comp.prop;

        // Bound constraint
        if prop.affine && 
           prop.a.len() == 1 && 
           *prop.a.values().next().unwrap() == 1. && 
           *self.kind() != ConstraintKind::Equal {
            match self.kind() {
                ConstraintKind::LessEqual => {    // x + b <= 0
                    u.push((prop.a.keys().next().unwrap().clone(),
                            -prop.b,
                            self.clone()))
                },
                ConstraintKind::GreaterEqual => { // x + b >= 0
                    l.push((prop.a.keys().next().unwrap().clone(),
                            -prop.b,
                            self.clone()))
                },
                _ => panic!("unexpected constraint type"),
            }
        }

        // Affine constraint
        else if prop.affine {

            // a^Tx + b == 0
            if *self.kind() == ConstraintKind::Equal {
                for (x, val) in prop.a.iter() {
                    a.push((*arow, x.clone(), *val)); 
                }
                b.push(-prop.b);
                ca.push(self.clone());
                *arow += 1;
            }

            // a^Tx + b - s == 0 and s <= 0 or s >= 0
            else {
                let s = VariableScalar::new_continuous(format!("_s_a{}_", *arow).as_str());
                for (x, val) in prop.a.iter() {
                    a.push((*arow, x.clone(), *val)); 
                }
                a.push((*arow, s.clone(), -1.)); 
                b.push(-prop.b);
                ca.push(self.clone());
                match self.kind() {
                    ConstraintKind::LessEqual => {
                        u.push((s.clone(), 0., self.clone()))
                    },
                    ConstraintKind::GreaterEqual => {
                        l.push((s.clone(), 0., self.clone()))
                    },
                    _ => panic!("unexpected constraint type"),
                }
                prop.a.insert(s.clone(), -1.);
                *arow += 1;
            }
        }

        // Nonlinear constraint
        else {

            // H
            let mut hh: Vec<(Node, Node, Node)> = Vec::new();
            for (v1, v2, e) in comp.hphi.iter() {
                hh.push((v1.clone(), v2.clone(), e.clone()))
            }
            h.push(hh);   

            // f(x) = 0
            if *self.kind() == ConstraintKind::Equal {
                f.push(comp.phi.clone());
                cj.push(self.clone());
                for (x, e) in comp.gphi.iter() {
                    j.push((*jrow, x.clone(), e.clone()));
                }
                *jrow += 1;
            }

            // f(x) - s == 0 and s <= 0 or s >= 0
            else {
                let s = VariableScalar::new_continuous(format!("_s_j{}_", *jrow).as_str());
                f.push(comp.phi-&s);
                cj.push(self.clone());
                for (x, e) in comp.gphi.iter() {
                    j.push((*jrow, x.clone(), e.clone()));
                }
                j.push((*jrow, s.clone(), ConstantScalar::new(-1.)));
                match self.kind() {
                    ConstraintKind::LessEqual => {
                        u.push((s.clone(), 0., self.clone()))
                    },
                    ConstraintKind::GreaterEqual => {
                        l.push((s.clone(), 0., self.clone()))
                    },
                    _ => panic!("unexpected constraint type"),
                }
                prop.a.insert(s.clone(), -1.);
                *jrow += 1;
            }
        }

        // Return
        ConstraintStdComp {
            ca: ca,
            cj: cj,
            a: a,
            b: b,
            f: f,
            j: j,
            h: h,
            u: u,
            l: l,
            prop: vec![prop],
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::model::node_cmp::NodeCmp;
    use crate::model::variable::VariableScalar;

    #[test]
    fn constr_std_comp_u_bound() {

        let x = VariableScalar::new_continuous("x");

        let c1 = (&x).leq(3.);
        let mut arow: usize = 1;
        let mut jrow: usize = 2;
        let comp1 = c1.std_components(&mut arow, &mut jrow);

        assert_eq!(comp1.ca.len(), 0);
        assert_eq!(comp1.cj.len(), 0);
        assert_eq!(comp1.a.len(), 0);
        assert_eq!(comp1.b.len(), 0);
        assert_eq!(comp1.f.len(), 0);
        assert_eq!(comp1.j.len(), 0);
        assert_eq!(comp1.h.len(), 0);
        assert_eq!(comp1.u.len(), 1);
        let (cx, cval, c) = &comp1.u[0];
        assert_eq!(*cx, x);
        assert_eq!(*cval, 3.);
        assert_eq!(*c, c1);
        assert_eq!(comp1.l.len(), 0);
        assert_eq!(arow, 1);
        assert_eq!(jrow, 2);
    }

    #[test]
    fn constr_std_comp_l_bound() {

        let x = VariableScalar::new_continuous("x");

        let c1 = (&x).geq(-4.);
        let mut arow: usize = 1;
        let mut jrow: usize = 2;
        let comp1 = c1.std_components(&mut arow, &mut jrow);

        assert_eq!(comp1.ca.len(), 0);
        assert_eq!(comp1.cj.len(), 0);
        assert_eq!(comp1.a.len(), 0);
        assert_eq!(comp1.b.len(), 0);
        assert_eq!(comp1.f.len(), 0);
        assert_eq!(comp1.j.len(), 0);
        assert_eq!(comp1.h.len(), 0);
        assert_eq!(comp1.u.len(), 0);
        assert_eq!(comp1.l.len(), 1);
        let (cx, cval, c) = &comp1.l[0];
        assert_eq!(*cx, x);
        assert_eq!(*cval, -4.);
        assert_eq!(*c, c1);
        assert_eq!(arow, 1);
        assert_eq!(jrow, 2);
    }

    #[test]
    fn constr_std_comp_affine_eq() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let c1 = (3.*&x + 4.*&y + 6.).equal(5.);
        let mut arow: usize = 1;
        let mut jrow: usize = 2;
        let comp1 = c1.std_components(&mut arow, &mut jrow);

        assert_eq!(comp1.ca.len(), 1);
        assert_eq!(comp1.ca[0], c1);
        assert_eq!(comp1.cj.len(), 0);
        assert_eq!(comp1.a.len(), 2);
        assert_eq!(comp1.b.len(), 1);
        let mut counter = 0_usize;
        for (row, col, val) in comp1.a.iter() {
            if *col == x {
                assert_eq!(*row, 1);
                assert_eq!(*val, 3.);
                counter += 1;
            }
            else if *col == y {
                assert_eq!(*row, 1);
                assert_eq!(*val, 4.);
                counter += 1;
            }
            else {
                panic!("unexpected variable");
            }
        }
        assert_eq!(counter, 2);
        assert_eq!(comp1.b[0], -1.);
        assert_eq!(comp1.f.len(), 0);
        assert_eq!(comp1.j.len(), 0);
        assert_eq!(comp1.h.len(), 0);
        assert_eq!(comp1.u.len(), 0);
        assert_eq!(comp1.l.len(), 0);
        assert_eq!(arow, 2);
        assert_eq!(jrow, 2);
    }

    #[test]
    fn constr_std_comp_affine_leq() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let c1 = (3.*&x + 4.*&y + 6.).leq(5.);
        let mut arow: usize = 1;
        let mut jrow: usize = 2;
        let comp1 = c1.std_components(&mut arow, &mut jrow);

        assert_eq!(comp1.ca.len(), 1);
        assert_eq!(comp1.ca[0], c1);
        assert_eq!(comp1.cj.len(), 0);
        assert_eq!(comp1.a.len(), 3);
        assert_eq!(comp1.b.len(), 1);
        let mut counter = 0_usize;
        for (row, col, val) in comp1.a.iter() {
            if *col == x {
                assert_eq!(*row, 1);
                assert_eq!(*val, 3.);
                counter += 1;
            }
            else if *col == y {
                assert_eq!(*row, 1);
                assert_eq!(*val, 4.);
                counter += 1;
            }
            else if (*col).name() == "_s_a1_" {
                assert_eq!(*row, 1);
                assert_eq!(*val, -1.);
                counter += 1;
            }
            else {
                panic!("unexpected variable");
            }
        }
        assert_eq!(counter, 3);
        assert_eq!(comp1.b[0], -1.);
        assert_eq!(comp1.f.len(), 0);
        assert_eq!(comp1.j.len(), 0);
        assert_eq!(comp1.h.len(), 0);
        assert_eq!(comp1.u.len(), 1);
        let (var, val, c) = &comp1.u[0];
        assert_eq!((*var).name(), "_s_a1_");
        assert_eq!(*val, 0.);
        assert_eq!(*c, c1);
        assert_eq!(comp1.l.len(), 0);
        assert_eq!(arow, 2);
        assert_eq!(jrow, 2);
    }

    #[test]
    fn constr_std_comp_affine_geq() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let c1 = (3.*&x + 4.*&y + 6.).geq(5.);
        let mut arow: usize = 1;
        let mut jrow: usize = 2;
        let comp1 = c1.std_components(&mut arow, &mut jrow);

        assert_eq!(comp1.ca.len(), 1);
        assert_eq!(comp1.ca[0], c1);
        assert_eq!(comp1.cj.len(), 0);
        assert_eq!(comp1.a.len(), 3);
        assert_eq!(comp1.b.len(), 1);
        let mut counter = 0_usize;
        for (row, col, val) in comp1.a.iter() {
            if *col == x {
                assert_eq!(*row, 1);
                assert_eq!(*val, 3.);
                counter += 1;
            }
            else if *col == y {
                assert_eq!(*row, 1);
                assert_eq!(*val, 4.);
                counter += 1;
            }
            else if (*col).name() == "_s_a1_" {
                assert_eq!(*row, 1);
                assert_eq!(*val, -1.);
                counter += 1;
            }
            else {
                panic!("unexpected variable");
            }
        }
        assert_eq!(counter, 3);    
        assert_eq!(comp1.b[0], -1.);
        assert_eq!(comp1.f.len(), 0);
        assert_eq!(comp1.j.len(), 0);
        assert_eq!(comp1.h.len(), 0);
        assert_eq!(comp1.u.len(), 0);
        assert_eq!(comp1.l.len(), 1);
        let (var, val, c) = &comp1.l[0];
        assert_eq!((*var).name(), "_s_a1_");
        assert_eq!(*val, 0.);
        assert_eq!(*c, c1);
        assert_eq!(arow, 2);
        assert_eq!(jrow, 2);
    }

    #[test]
    fn constr_std_comp_nonlinear_eq() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let c1 = (3.*&x*&x + 4.*&x*&y + 7.*&y*&y + 8.).equal(5.);
        let mut arow: usize = 1;
        let mut jrow: usize = 2;
        let comp1 = c1.std_components(&mut arow, &mut jrow);

        assert_eq!(comp1.ca.len(), 0);
        assert_eq!(comp1.cj.len(), 1);
        assert_eq!(comp1.a.len(), 0);
        assert_eq!(comp1.b.len(), 0);
        assert_eq!(comp1.f.len(), 1);
        assert_eq!(format!("{}", comp1.f[0]),
                   "3*x*x + 4*x*y + 7*y*y + 3");
        assert_eq!(comp1.j.len(), 2);
        let mut counter: usize = 0;
        for (row, col, val) in comp1.j.iter() {
            if *col == x {
                assert_eq!(*row, 2);
                assert_eq!(format!("{}", *val),
                           "y*4 + 3*x + x*3");
                counter += 1;
            }
            else if *col == y {
                assert_eq!(*row, 2);
                assert_eq!(format!("{}", *val),
                           "7*y + y*7 + 4*x");
                counter += 1;
            }
            else {
                panic!("unexpected variable");
            }
        }
        assert_eq!(counter, 2);
        assert_eq!(comp1.h.len(), 1);
        assert_eq!(comp1.h[0].len(), 3);
        let mut counter: usize = 0;
        for (row, col, val) in comp1.h[0].iter() {
            if *row == x && *col == x {
                assert!((*val).is_constant_with_value(6.));
                counter += 1;
            }
            else if *row == x && *col == y {
                assert!((*val).is_constant_with_value(4.));
                counter += 1;
            }
            else if *row == y && *col == x{
                assert!((*val).is_constant_with_value(4.));
                counter += 1;
            }
            else if *row == y && *col == y {
                assert!((*val).is_constant_with_value(14.));
                counter += 1;
            }
            else {
                panic!("unexpected variable pair")
            }
        }
        assert_eq!(counter, 3);
        assert_eq!(comp1.u.len(), 0);
        assert_eq!(comp1.l.len(), 0);
        assert_eq!(arow, 1);
        assert_eq!(jrow, 3);
    }

    #[test]
    fn constr_std_comp_nonlinear_leq() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let c1 = (3.*&x*&x + 4.*&x*&y + 7.*&y*&y + 8.).leq(5.);
        let mut arow: usize = 1;
        let mut jrow: usize = 2;
        let comp1 = c1.std_components(&mut arow, &mut jrow);

        assert_eq!(comp1.ca.len(), 0);
        assert_eq!(comp1.cj.len(), 1);
        assert_eq!(comp1.a.len(), 0);
        assert_eq!(comp1.b.len(), 0);
        assert_eq!(comp1.f.len(), 1);
        assert_eq!(format!("{}", comp1.f[0]),
                   "3*x*x + 4*x*y + 7*y*y + -1*_s_j2_ + 3");
        assert_eq!(comp1.j.len(), 3);
        let mut counter: usize = 0;
        for (row, col, val) in comp1.j.iter() {
            if *col == x {
                assert_eq!(*row, 2);
                assert_eq!(format!("{}", *val),
                           "y*4 + 3*x + x*3");
                counter += 1;
            }
            else if *col == y {
                assert_eq!(*row, 2);
                assert_eq!(format!("{}", *val),
                           "7*y + y*7 + 4*x");
                counter += 1;
            }
            else if (*col).name() == "_s_j2_"  {
                assert_eq!(*row, 2);
                assert!((*val).is_constant_with_value(-1.));
                counter += 1;
            }
            else {
                panic!("unexpected variable");
            }
        }
        assert_eq!(counter, 3);
        assert_eq!(comp1.h.len(), 1);
        assert_eq!(comp1.h[0].len(), 3);
        let mut counter: usize = 0;
        for (row, col, val) in comp1.h[0].iter() {
            if *row == x && *col == x {
                assert!((*val).is_constant_with_value(6.));
                counter += 1;
            }
            else if *row == x && *col == y {
                assert!((*val).is_constant_with_value(4.));
                counter += 1;
            }
            else if *row == y && *col == x{
                assert!((*val).is_constant_with_value(4.));
                counter += 1;
            }
            else if *row == y && *col == y {
                assert!((*val).is_constant_with_value(14.));
                counter += 1;
            }
            else {
                panic!("unexpected variable pair")
            }
        }
        assert_eq!(counter, 3);
        assert_eq!(comp1.u.len(), 1);
        let (var, val, c) = comp1.u.iter().next().unwrap();
        assert_eq!((*var).name(), "_s_j2_");
        assert_eq!(*val, 0.);
        assert_eq!(*c, c1);
        assert_eq!(comp1.l.len(), 0);
        assert_eq!(arow, 1);
        assert_eq!(jrow, 3);
    }

    #[test]
    fn constr_std_comp_nonlinear_geq() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let c1 = (3.*&x*&x + 4.*&x*&y + 7.*&y*&y + 8.).geq(5.);
        let mut arow: usize = 1;
        let mut jrow: usize = 2;
        let comp1 = c1.std_components(&mut arow, &mut jrow);

        assert_eq!(comp1.ca.len(), 0);
        assert_eq!(comp1.cj.len(), 1);
        assert_eq!(comp1.a.len(), 0);
        assert_eq!(comp1.b.len(), 0);
        assert_eq!(comp1.f.len(), 1);
        assert_eq!(format!("{}", comp1.f[0]),
                   "3*x*x + 4*x*y + 7*y*y + -1*_s_j2_ + 3");
        assert_eq!(comp1.j.len(), 3);
        let mut counter: usize = 0;
        for (row, col, val) in comp1.j.iter() {
            if *col == x {
                assert_eq!(*row, 2);
                assert_eq!(format!("{}", *val),
                           "y*4 + 3*x + x*3");
                counter += 1;
            }
            else if *col == y {
                assert_eq!(*row, 2);
                assert_eq!(format!("{}", *val),
                           "7*y + y*7 + 4*x");
                counter += 1;
            }
            else if (*col).name() == "_s_j2_" {
                assert_eq!(*row, 2);
                assert!((*val).is_constant_with_value(-1.));
                counter += 1;
            }
            else {
                panic!("unexpected variable");
            }
        }
        assert_eq!(counter, 3);
        assert_eq!(comp1.h.len(), 1);
        assert_eq!(comp1.h[0].len(), 3);
        let mut counter: usize = 0;
        for (row, col, val) in comp1.h[0].iter() {
            if *row == x && *col == x {
                assert!((*val).is_constant_with_value(6.));
                counter += 1;
            }
            else if *row == x && *col == y {
                assert!((*val).is_constant_with_value(4.));
                counter += 1;
            }
            else if *row == y && *col == x{
                assert!((*val).is_constant_with_value(4.));
                counter += 1;
            }
            else if *row == y && *col == y {
                assert!((*val).is_constant_with_value(14.));
                counter += 1;
            }
            else {
                panic!("unexpected variable pair")
            }
        }
        assert_eq!(counter, 3);
        assert_eq!(comp1.u.len(), 0);
        assert_eq!(comp1.l.len(), 1);
        let (var, val, c) = comp1.l.iter().next().unwrap();
        assert_eq!((*var).name(), "_s_j2_");
        assert_eq!(*val, 0.);
        assert_eq!(*c, c1);
        assert_eq!(arow, 1);
        assert_eq!(jrow, 3);
    }
}