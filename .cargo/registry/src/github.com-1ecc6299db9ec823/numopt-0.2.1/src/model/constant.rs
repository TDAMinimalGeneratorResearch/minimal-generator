//! Constant expression.

use std::fmt;
use std::rc::Rc;

use super::node::Node;
use super::node_base::NodeBase;

/// Constant scalar.
pub struct ConstantScalar {
    value: f64,
}

impl ConstantScalar {

    /// Creates new scalar constant expression node.
    pub fn new(value: f64) -> Node {
        Node::ConstantScalar(Rc::new(
            Self {
                value: value,
            }
        ))
    }

    /// Gets value of constant.
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Creates new scalar constant with value zero.
    pub fn zero() -> Node {
        ConstantScalar::new(0.)
    }
}

impl NodeBase for ConstantScalar {

    fn partial(&self, _arg: &Node) -> Node { ConstantScalar::new(0.) }
}

impl fmt::Display for ConstantScalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

#[cfg(test)]
mod tests {

    use crate::model::node_base::NodeBase;
    use crate::model::node_std::NodeStd;
    use crate::model::node_diff::NodeDiff;
    use crate::model::variable::VariableScalar;
    use crate::model::constant::ConstantScalar;

    #[test]
    fn constant_zero() {

        let z = ConstantScalar::zero();
        assert!(z.is_constant_with_value(0.));
    }

    #[test]
    fn constant_partial() {

        let x = VariableScalar::new_continuous("x");
        let c = ConstantScalar::new(3.);

        let z1 = c.partial(&x);
        assert!(z1.is_constant_with_value(0.));

        let z2 = c.partial(&c);
        assert!(z2.is_constant_with_value(0.));
    }

    #[test]
    fn constant_derivative() {

        let x = VariableScalar::new_continuous("x");
        let c = ConstantScalar::new(3.);

        let z1 = c.derivative(&x);
        assert!(z1.is_constant_with_value(0.));
    }

    #[test]
    fn constant_std_properties() {

        let c = ConstantScalar::new(5.);
        let p = c.std_properties();
        assert!(p.affine);
        assert_eq!(p.b, 5.);
        assert!(p.a.is_empty());
    }
}