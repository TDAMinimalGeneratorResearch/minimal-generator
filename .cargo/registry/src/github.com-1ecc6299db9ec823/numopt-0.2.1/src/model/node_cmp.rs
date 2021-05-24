//! Trait for comparing expression nodes and constructing 
//! optimization constraints. 

use num_traits::cast::ToPrimitive;

use crate::model::node::Node;
use crate::model::constant::ConstantScalar;
use crate::model::constraint::{Constraint, ConstraintKind};

/// Trait for comparing expression nodes.
pub trait NodeCmp<T> {

    /// Creates equality constraint and tags it.
    fn equal_and_tag(&self, other: T, tag: &str) -> Constraint;

    /// Creates equality constraint.
    fn equal(&self, other: T) -> Constraint { self.equal_and_tag(other, "") }

    /// Creates greater-than-or-equal constraint and tags it.
    fn geq_and_tag(&self, other: T, tag: &str) -> Constraint;

    /// Creates greater-than-or-equal constraint.
    fn geq(&self, other: T) -> Constraint { self.geq_and_tag(other, "") }

    /// Creates less-than-or-equal constraint and tags it.
    fn leq_and_tag(&self, other: T, tag: &str) -> Constraint;

    /// Creates less-than-or-equal constraint.
    fn leq(&self, other: T) -> Constraint { self.leq_and_tag(other, "") }
}

macro_rules! impl_node_cmp_scalar {
    ($x: ty, $y: ty) => {
        impl NodeCmp<$y> for $x {
    
            fn equal_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(self.clone(),
                                ConstraintKind::Equal,
                                ConstantScalar::new(other.to_f64().unwrap()),
                                tag)
            }

            fn geq_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(self.clone(),
                                ConstraintKind::GreaterEqual,
                                ConstantScalar::new(other.to_f64().unwrap()),
                                tag)
            }

            fn leq_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(self.clone(),
                                ConstraintKind::LessEqual,
                                ConstantScalar::new(other.to_f64().unwrap()),
                                tag)
            }
        }
    };
}

impl_node_cmp_scalar!(Node, f64);

macro_rules! impl_node_cmp_node {
    ($x: ty, $y: ty) => {
        impl NodeCmp<$y> for $x {
    
            fn equal_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(self.clone(),
                                ConstraintKind::Equal,
                                other.clone(),
                                tag)
            }

            fn geq_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(self.clone(),
                                ConstraintKind::GreaterEqual,
                                other.clone(),
                                tag)
            }

            fn leq_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(self.clone(),
                                ConstraintKind::LessEqual,
                                other.clone(),
                                tag)
            }
        }
    };
}

impl_node_cmp_node!(Node, Node);
impl_node_cmp_node!(Node, &Node);

macro_rules! impl_scalar_cmp_node {
    ($x: ty, $y: ty) => {
        impl NodeCmp<$y> for $x {
    
            fn equal_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(ConstantScalar::new(self.to_f64().unwrap()),
                                ConstraintKind::Equal,
                                other.clone(),
                                tag)
            }

            fn geq_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(ConstantScalar::new(self.to_f64().unwrap()),
                                ConstraintKind::GreaterEqual,
                                other.clone(),
                                tag)
            }

            fn leq_and_tag(&self, other: $y, tag: &str) -> Constraint {
                Constraint::new(ConstantScalar::new(self.to_f64().unwrap()),
                                ConstraintKind::LessEqual,
                                other.clone(),
                                tag)
            }
        }
    };
}

impl_scalar_cmp_node!(f64, Node);
impl_scalar_cmp_node!(f64, &Node);

#[cfg(test)]
mod tests {

    use crate::model::node_cmp::NodeCmp;
    use crate::model::variable::VariableScalar;
    use crate::model::constant::ConstantScalar;

    #[test]
    fn node_cmp_node() {

        let x = VariableScalar::new_continuous("x");
        let c = ConstantScalar::new(5.);

        let z1 = x.equal(&c);
        assert_eq!(format!("{}", z1), "x == 5");

        let z2 = &x.leq(&c);
        assert_eq!(format!("{}", z2), "x <= 5");

        let z3 = x.geq(&x + 3.);
        assert_eq!(format!("{}", z3), "x >= x + 3");

        let z4 = &x.leq(5.*&x);
        assert_eq!(format!("{}", z4), "x <= 5*x");
    }

    #[test]
    fn node_cmp_scalar() {

        let x = VariableScalar::new_continuous("x");

        let z1 = x.equal(6.);
        assert_eq!(format!("{}", z1), "x == 6");

        let z2 = &x.equal(10.);
        assert_eq!(format!("{}", z2), "x == 10");

        let z3 = (&x + 11.).equal(12.);
        assert_eq!(format!("{}", z3), "x + 11 == 12");
    }

    #[test]
    fn scalar_cmp_node() {

        let x = VariableScalar::new_continuous("x");

        let z1 = 4_f64.equal(&x);
        assert_eq!(format!("{}", z1), "4 == x");

        let z2 = 4_f64.leq(&x + 3.);
        assert_eq!(format!("{}", z2), "4 <= x + 3");

        let z3 = 5_f64.geq(&x*5.);
        assert_eq!(format!("{}", z3), "5 >= x*5");
    }
}