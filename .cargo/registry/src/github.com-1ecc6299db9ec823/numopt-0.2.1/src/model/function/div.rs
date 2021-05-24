//! Divide function.

use std::fmt;
use std::rc::Rc;
use std::collections::HashMap;

use crate::model::node::Node;
use crate::model::node_base::NodeBase;
use crate::model::node_std::{NodeStd, NodeStdProp};
use crate::model::constant::ConstantScalar;

/// Divide function.
pub struct FunctionDiv {
    args: (Node, Node),
}

impl FunctionDiv {

    /// Creates new divide expression node.
    pub fn new(arg1: Node, arg2: Node) -> Node {
        Node::FunctionDiv(Rc::new(
            Self {
                args: (arg1, arg2),
            }
        ))
    }
}

impl NodeBase for FunctionDiv {

    fn arguments(&self) -> Vec<&Node> {
        vec![&self.args.0, &self.args.1]
    }

    fn partial(&self, arg: &Node) -> Node { 
        if self.args.0 == *arg {
            return 1./&self.args.1;
        }
        else if self.args.1 == *arg {
            return -&self.args.0/(&self.args.1*&self.args.1);
        }
        else {
            return ConstantScalar::new(0.);
        }
    }

    fn evaluate(&self, var_values: &HashMap<&Node, f64>) -> f64 { 
        self.args.0.evaluate(var_values)/self.args.1.evaluate(var_values)
    }
}

impl NodeStd for FunctionDiv {

    fn std_properties(&self) -> NodeStdProp {
        
        let p0 = self.args.0.std_properties();
        let p1 = self.args.1.std_properties();
        let affine = p0.affine && p1.a.is_empty();
        let b = p0.b/p1.b;
        let mut a: HashMap<Node, f64> = HashMap::new();
        for (key, val) in p0.a.iter() {
            a.insert(key.clone(), (*val)/p1.b);
        }
        for (key, _val) in p1.a.iter() {
            a.insert(key.clone(), 0.);
        }
        NodeStdProp {
            affine: affine,
            a: a,
            b: b,
        }
    }
}

impl<'a> fmt::Display for FunctionDiv {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s0 = match &self.args.0 {
            Node::FunctionAdd(x) => format!("({})", x),
            Node::FunctionDiv(x) => format!("({})", x),
            _ => format!("{}", self.args.0)
        };
        let s1 = match &self.args.1 {
            Node::FunctionAdd(x) => format!("({})", x),
            Node::FunctionMul(x) => format!("({})", x),
            Node::FunctionDiv(x) => format!("({})", x),
            _ => format!("{}", self.args.1)
        };
        write!(f, "{}/{}", s0, s1)
    }
}

#[cfg(test)]
mod tests {

    use crate::model::node_base::NodeBase;
    use crate::model::node_std::NodeStd;
    use crate::model::node_diff::NodeDiff;
    use crate::model::variable::VariableScalar;

    #[test]
    fn div_partial() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let w = VariableScalar::new_continuous("w");

        let z = &x/&y;

        let z1 = z.partial(&x);
        assert_eq!(format!("{}", z1), "1/y");

        let z2 = z.partial(&y);
        assert_eq!(format!("{}", z2), "-1*x/(y*y)");

        let z3 = z.partial(&w);
        assert!(z3.is_constant_with_value(0.));
    }

    #[test]
    fn div_derivative() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let z1 = &x/3.;
        let z1x = z1.derivative(&x);
        let z1y = z1.derivative(&y);
        assert!(z1x.is_constant_with_value(1./3.));
        assert!(z1y.is_constant_with_value(0.));

        let z2 = 4./&x;
        let z2x = z2.derivative(&x);
        let z2y = z2.derivative(&y);
        assert_eq!(format!("{}", z2x), "-4/(x*x)");
        assert!(z2y.is_constant_with_value(0.));

        let z3 = 5./-&y;
        let z3x = z3.derivative(&x);
        let z3y = z3.derivative(&y);
        assert!(z3x.is_constant_with_value(0.));
        assert_eq!(format!("{}", z3y), "(-5/(-1*y*-1*y))*-1");

        let z4 = 3.*&x/(&y - &x);
        let z4x = z4.derivative(&x);
        let z4y = z4.derivative(&y);
        assert_eq!(format!("{}", z4x), 
                  "(-1*3*x/((y + -1*x)*(y + -1*x)))*-1 + (1/(y + -1*x))*3"); 
        assert_eq!(format!("{}", z4y), 
                  "-1*3*x/((y + -1*x)*(y + -1*x))");

        let f1 = &x - 2.;
        let z5 = &f1/(&f1 + 3.);
        let z5x = z5.derivative(&x);
        assert_eq!(format!("{}", z5x),
                   "(-1*x + 2)/((x + 1)*(x + 1)) + 1/(x + 1)");
    }

    #[test]
    fn div_std_properties() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let z1 = 3./&x;
        let p1 = z1.std_properties();
        assert!(!p1.affine);
        assert_eq!(p1.a.len(), 1);
        assert!(p1.a.contains_key(&x));
        
        let z2 = (&x + 10.)/4.;
        let p2 = z2.std_properties();
        assert!(p2.affine);
        assert_eq!(p2.b, 10./4.);
        assert_eq!(p2.a.len(), 1);
        assert_eq!(*p2.a.get(&x).unwrap(), 1./4.);

        let z3 = 3./(&x + &y);
        let p3 = z3.std_properties();
        assert!(!p3.affine);
        assert_eq!(p3.a.len(), 2);
        assert!(p3.a.contains_key(&x));
        assert!(p3.a.contains_key(&y));

        let z4 = (4.*&x + 5. + &y)/10.;
        let p4 = z4.std_properties();
        assert!(p4.affine);
        assert_eq!(p4.b, 0.5);
        assert_eq!(p4.a.len(), 2);
        assert_eq!(*p4.a.get(&x).unwrap(), 4./10.);
        assert_eq!(*p4.a.get(&y).unwrap(), 1./10.);

        let z5 = (4.*&y)/(5. + &x);
        let p5 = z5.std_properties();
        assert!(!p5.affine);
        assert_eq!(p5.a.len(), 2);
        assert!(p5.a.contains_key(&x));
        assert!(p5.a.contains_key(&y));
    }
}


