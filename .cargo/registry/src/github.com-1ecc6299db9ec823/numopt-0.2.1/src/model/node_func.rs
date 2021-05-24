//! Trait that provides various functions of expression nodes.

use crate::model::node::Node;
use crate::model::constant::ConstantScalar;
use crate::model::function::cos::FunctionCos;
use crate::model::function::sin::FunctionSin;

/// Trait that provides various functions for expression nodes.
pub trait NodeFunc {

    /// Constructs cosine expression node.
    fn cos(&self) -> Node;

    /// Constructs sine expression node.
    fn sin(&self) -> Node;
}

impl NodeFunc for Node {

    fn cos(&self) -> Node {
        match self {
            Node::ConstantScalar(x) => {
                ConstantScalar::new(x.value().cos())
            },
            _ =>  FunctionCos::new(self.clone())  
        }
    }

    fn sin(&self) -> Node {
        match self {
            Node::ConstantScalar(x) => {
                ConstantScalar::new(x.value().sin())
            },
            _ =>  FunctionSin::new(self.clone())  
        }
    }
}

#[cfg(test)]
mod tests {

    use maplit::hashmap;

    use super::*;
    use crate::model::node_base::NodeBase;
    use crate::model::variable::VariableScalar;
    use crate::model::constant::ConstantScalar;

    #[test]
    fn node_cos() {

        let x = VariableScalar::new_continuous("x");
        let c = ConstantScalar::new(5.);

        let var_values = hashmap!{ &x => 3. }; 

        let z1 = x.cos();
        assert_eq!(format!("{}", z1), "cos(x)");
        assert_eq!(z1.evaluate(&var_values), 3_f64.cos());

        let z2 = (3.*&x + 5.).cos();
        assert_eq!(format!("{}", z2), "cos(3*x + 5)");
        assert_eq!(z2.evaluate(&var_values), (3.*3. + 5_f64).cos());

        let z3 = c.cos();
        assert!(z3.is_constant_with_value(5_f64.cos()));
    }

    #[test]
    fn node_sin() {

        let x = VariableScalar::new_continuous("x");
        let c = ConstantScalar::new(5.);

        let var_values = hashmap!{ &x => 3. }; 

        let z1 = x.sin();
        assert_eq!(format!("{}", z1), "sin(x)");
        assert_eq!(z1.evaluate(&var_values), 3_f64.sin());

        let z2 = (3.*&x + 5.).sin();
        assert_eq!(format!("{}", z2), "sin(3*x + 5)");
        assert_eq!(z2.evaluate(&var_values), (3.*3. + 5_f64).sin());

        let z3 = c.sin();
        assert!(z3.is_constant_with_value(5_f64.sin()));
    }
}