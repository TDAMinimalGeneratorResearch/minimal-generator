//! Optimization variables.

use std::fmt;
use std::rc::Rc;

use super::node::Node;
use super::node_base::NodeBase;
use super::constant::ConstantScalar;

/// Optimization variable kind.
pub enum VariableKind {
    VarContinuous,
    VarInteger,
}

/// Scalar optimization variable.
pub struct VariableScalar {
    name: String,
    kind: VariableKind,
}

impl VariableScalar {

    /// Determines whether optimization variable is a continuous variable.
    pub fn is_continuous(&self) -> bool {
        match self.kind {
            VariableKind::VarContinuous => true,
            _ => false,
        }
    }

    /// Determines whether optimization variable is an integer variable.
    pub fn is_integer(&self) -> bool {
        match self.kind {
            VariableKind::VarInteger => true,
            _ => false,
        }
    }

    /// Gets optimization variable name.
    pub fn name(&self) -> &str { self.name.as_ref() }

    /// Creates a new optimization variable expression node.
    pub fn new(name: &str, kind: VariableKind) -> Node {
        Node::VariableScalar(Rc::new(
            Self {
                name: name.to_string(),
                kind: kind,
            }
        ))
    }

    /// Creates a new continuous optimization variable expression node.
    pub fn new_continuous(name: &str) -> Node {
        VariableScalar::new(name, VariableKind::VarContinuous)
    }

    /// Creates a new integer optimization variable expression node.
    pub fn new_integer(name: &str) -> Node {
        VariableScalar::new(name, VariableKind::VarInteger)
    }
}

impl NodeBase for VariableScalar {

    fn partial(&self, arg: &Node) -> Node { 
        match arg {
            Node::VariableScalar(x) => {
                if self as *const VariableScalar == x.as_ref() {
                    ConstantScalar::new(1.)       
                }
                else {
                    ConstantScalar::new(0.)       
                }
            }
            _ => ConstantScalar::new(0.)  
        }
    }
}

impl<'a> fmt::Display for VariableScalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name)
    }
}

#[cfg(test)]
mod tests {

    use crate::model::node::Node;
    use crate::model::node_base::NodeBase;
    use crate::model::node_std::NodeStd;
    use crate::model::node_diff::NodeDiff;
    use crate::model::variable::VariableScalar;

    #[test]
    fn var_construction() {

        let x = VariableScalar::new_continuous("x");
        assert_eq!(x.name(), "x");
        match x {
            Node::VariableScalar(xx) => {
                assert!(xx.is_continuous());
                assert!(!xx.is_integer());
            },
            _ => panic!("construction failed"),
        }

        let y = VariableScalar::new_integer("y");
        assert_eq!(y.name(), "y");
        match y {
            Node::VariableScalar(yy) => {
                assert!(yy.is_integer());
                assert!(!yy.is_continuous());
            },
            _ => panic!("construction failed"),
        }
    }

    #[test]
    fn var_partial() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let z1 = x.partial(&x);
        assert!(z1.is_constant_with_value(1.));

        let z2 = x.partial(&y);
        assert!(z2.is_constant_with_value(0.));
    }

    #[test]
    fn var_derivative() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let z1 = x.derivative(&y);
        assert!(z1.is_constant_with_value(0.));

        let z2 = x.derivative(&x);
        assert!(z2.is_constant_with_value(1.));
    }

    #[test]
    fn var_std_properties() {

        let x = VariableScalar::new_integer("x");
        let p = x.std_properties();
        assert!(p.affine);
        assert_eq!(p.b, 0.);
        assert_eq!(p.a.len(), 1);
        assert_eq!(*p.a.get(&x).unwrap(), 1.);
    }
}