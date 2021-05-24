//! Multiply function.

use std::fmt;
use std::rc::Rc;
use std::collections::HashMap;

use crate::model::node::Node;
use crate::model::node_base::NodeBase;
use crate::model::node_std::{NodeStd, NodeStdProp};
use crate::model::constant::ConstantScalar;

/// Multiply function.
pub struct FunctionMul {
    args: (Node, Node),
}

impl FunctionMul {

    /// Creates new multiply expression node.
    pub fn new(arg1: Node, arg2: Node) -> Node {
        Node::FunctionMul(Rc::new(
            Self {
                args: (arg1, arg2),
            }
        ))
    }
}

impl NodeBase for FunctionMul {

    fn arguments(&self) -> Vec<&Node> {
        vec![&self.args.0, &self.args.1]
    }

    fn partial(&self, arg: &Node) -> Node { 
        if self.args.0 == *arg {
            return self.args.1.clone();
        }
        else if self.args.1 == *arg {
            return self.args.0.clone();
        }
        else {
            return ConstantScalar::new(0.);
        }
    }

    fn evaluate(&self, var_values: &HashMap<&Node, f64>) -> f64 { 
        self.args.0.evaluate(var_values)*self.args.1.evaluate(var_values)
    }
}

impl NodeStd for FunctionMul {

    fn std_properties(&self) -> NodeStdProp {
        
        let p0 = self.args.0.std_properties();
        let p1 = self.args.1.std_properties();
        let affine = (p0.affine && p1.a.is_empty()) || 
                     (p1.affine && p0.a.is_empty());
        let b = p0.b*p1.b;
        let mut a: HashMap<Node, f64> = HashMap::new();
        for (key, val) in p0.a.iter() {
            a.insert(key.clone(), (*val)*p1.b);
        }
        for (key, val) in p1.a.iter() {
            a.insert(key.clone(), (*val)*p0.b);
        }
        NodeStdProp {
            affine: affine,
            a: a,
            b: b,
        }
    }
}

impl<'a> fmt::Display for FunctionMul {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s0 = match &self.args.0 {
            Node::FunctionAdd(x) => format!("({})", x),
            Node::FunctionDiv(x) => format!("({})", x),
            _ => format!("{}", self.args.0)
        };
        let s1 = match &self.args.1 {
            Node::FunctionAdd(x) => format!("({})", x),
            _ => format!("{}", self.args.1)
        };
        write!(f, "{}*{}", s0, s1)
    }
}

#[cfg(test)]
mod tests {

    use crate::model::node_base::NodeBase;
    use crate::model::node_std::NodeStd;
    use crate::model::node_diff::NodeDiff;
    use crate::model::variable::VariableScalar;

    #[test]
    fn mul_partial() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let w = VariableScalar::new_continuous("w");

        let z = &x*&y;

        let z1 = z.partial(&x);
        assert_eq!(z1, y);

        let z2 = z.partial(&y);
        assert_eq!(z2, x);

        let z3 = z.partial(&w);
        assert!(z3.is_constant_with_value(0.));
    }

    #[test]
    fn mul_derivative() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let z1 = 3.*&x;
        let z1x = z1.derivative(&x);
        let z1y = z1.derivative(&y);
        assert!(z1x.is_constant_with_value(3.));
        assert!(z1y.is_constant_with_value(0.));

        let z2 = &x*&y;
        let z2x = z2.derivative(&x);
        let z2y = z2.derivative(&y);
        assert_eq!(z2x, y);
        assert_eq!(z2y, x);
        
        let z3 = &y*(&x - 3. - &y*&y);
        let z3x = z3.derivative(&x);
        let z3y = z3.derivative(&y);
        assert_eq!(z3x, y);
        assert_eq!(format!("{}", z3y), "y*-1*y + y*-1*y + x + -1*y*y + -3");

        let f1 = 3.*&x;
        let z4 = &f1*(&f1*&y);
        let z4x = z4.derivative(&x);
        let z4y = z4.derivative(&y);
        assert_eq!(format!("{}", z4x), "3*x*y*3 + 3*x*y*3");
        assert_eq!(format!("{}", z4y), "3*x*3*x");
    }

    #[test]
    fn mul_std_properties() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let z1 = 3.*&x;
        let p1 = z1.std_properties();
        assert!(p1.affine);
        assert_eq!(p1.b, 0.);
        assert_eq!(p1.a.len(), 1);
        assert_eq!(*p1.a.get(&x).unwrap(), 3.);
        
        let z2 = &x*4.;
        let p2 = z2.std_properties();
        assert!(p2.affine);
        assert_eq!(p2.b, 0.);
        assert_eq!(p2.a.len(), 1);
        assert_eq!(*p2.a.get(&x).unwrap(), 4.);

        let z3 = 3.*(&x + 3.);
        let p3 = z3.std_properties();
        assert!(p3.affine);
        assert_eq!(p3.b, 9.);
        assert_eq!(p3.a.len(), 1);
        assert_eq!(*p3.a.get(&x).unwrap(), 3.);

        let z4 = (4.*&x + 5.)*10.;
        let p4 = z4.std_properties();
        assert!(p4.affine);
        assert_eq!(p4.b, 50.);
        assert_eq!(p4.a.len(), 1);
        assert_eq!(*p4.a.get(&x).unwrap(), 40.);

        let z5 = (4.*&y)*(5. + &x);
        let p5 = z5.std_properties();
        assert!(!p5.affine);
        assert_eq!(p5.a.len(), 2);
        assert!(p5.a.contains_key(&x));
        assert!(p5.a.contains_key(&y));
    }
}

