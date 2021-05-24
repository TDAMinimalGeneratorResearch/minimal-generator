//! Add function.

use std::fmt;
use std::rc::Rc;
use std::collections::HashMap;

use crate::model::node::Node;
use crate::model::node_base::NodeBase;
use crate::model::node_std::{NodeStd, NodeStdProp};
use crate::model::constant::ConstantScalar;

/// Add function.
pub struct FunctionAdd {
    args: Vec<Node>,
}

impl FunctionAdd {

    /// Creates new add expression node.
    pub fn new(args: Vec<Node>) -> Node {

        assert!(args.len() >= 2);
        Node::FunctionAdd(Rc::new(
            Self {
                args: args,
            }
        ))
    }
}

impl NodeBase for FunctionAdd {

    fn arguments(&self) -> Vec<&Node> {
        self.args.iter().collect()
    }

    fn partial(&self, arg: &Node) -> Node { 
        for a in &self.args {
            if *a == *arg {
                return ConstantScalar::new(1.);
            } 
        }
        ConstantScalar::new(0.)
    }

    fn evaluate(&self, var_values: &HashMap<&Node, f64>) -> f64 { 
        self.args.iter().map(|x| x.evaluate(var_values)).sum()
    }
}

impl NodeStd for FunctionAdd {

    fn std_properties(&self) -> NodeStdProp {
        
        let mut affine = true;
        let mut a: HashMap<Node, f64> = HashMap::new();
        let mut b = 0_f64;
        for arg in self.arguments().iter() {
            let p = arg.std_properties();
            affine &= p.affine;
            b += p.b;
            for (key, val) in p.a.iter() {
                if let Some(x) = a.get_mut(key) {
                    *x += *val;
                }
                else {
                    a.insert(key.clone(), *val);
                }
            }
        }
        NodeStdProp {
            affine: affine,
            a: a,
            b: b,
        }
    }
}

impl<'a> fmt::Display for FunctionAdd {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = self.args.len();
        for i in 0..n {
            if i < n-1 {
                write!(f, "{} + ", self.args[i]).unwrap();
            }
            else {
                write!(f, "{}", self.args[i]).unwrap();
            }
        };
        Ok(())
    }
}

#[cfg(test)]
mod tests {

    use maplit::hashmap;

    use crate::model::node_base::NodeBase;
    use crate::model::node_std::NodeStd;
    use crate::model::node_func::NodeFunc;
    use crate::model::node_diff::NodeDiff;
    use crate::model::variable::VariableScalar;

    #[test]
    fn add_partial() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let w = VariableScalar::new_continuous("w");

        let z = &x + &y; 

        let z1 = z.partial(&x);
        assert!(z1.is_constant_with_value(1.));

        let z2 = z.partial(&y);
        assert!(z2.is_constant_with_value(1.));

        let z3 = z.partial(&w);
        assert!(z3.is_constant_with_value(0.));

        let zz = &x + 2.;
        let f = &y + &zz;

        let z4 = f.partial(&x);
        assert!(z4.is_constant_with_value(1.));

        let z5 = f.partial(&y);
        assert!(z5.is_constant_with_value(1.));

        let z6 = f.partial(&zz);
        assert!(z6.is_constant_with_value(0.));
    }

    #[test]
    fn add_derivative() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let z1 = &x + 1.;
        let z1x = z1.derivative(&x);
        let z1y = z1.derivative(&y);
        assert!(z1x.is_constant_with_value(1.));
        assert!(z1y.is_constant_with_value(0.));

        let z2 = &x + &y;
        let z2x = z2.derivative(&x);
        let z2y = z2.derivative(&y);
        assert!(z2x.is_constant_with_value(1.));
        assert!(z2y.is_constant_with_value(1.));

        let z3 = (&x + 1.) + (&x + 3.) + (&y + (&x + 5.));
        let z3x = z3.derivative(&x);
        let z3y = z3.derivative(&y);
        assert!(z3x.is_constant_with_value(3.));
        assert!(z3y.is_constant_with_value(1.));

        let z4 = &x + &x;
        let z4x = z4.derivative(&x);
        let z4y = z4.derivative(&y);
        assert!(z4x.is_constant_with_value(2.));
        assert!(z4y.is_constant_with_value(0.));

        let f1 = &x + 1. + &y;
        let z5 = &f1 + &f1;
        let z5x = z5.derivative(&x);
        let z5y = z5.derivative(&y);
        assert_eq!(z5.evaluate(&hashmap!{ &x => 3., &y => 4. }), 2.*(3.+1.+4.));
        assert!(z5x.is_constant_with_value(2.));
        assert!(z5y.is_constant_with_value(2.));
    }

    #[test]
    fn add_std_properties() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let z = VariableScalar::new_continuous("z");

        let z1 = &x + &x;
        let p1 = z1.std_properties();
        assert!(p1.affine);
        assert_eq!(p1.b, 0.);
        assert_eq!(p1.a.len(), 1);
        assert_eq!(*p1.a.get(&x).unwrap(), 2.);

        let z2 = &x + 8. + &y + 40. + &x;
        let p2 = z2.std_properties();
        assert!(p2.affine);
        assert_eq!(p2.b, 48.);
        assert_eq!(p2.a.len(), 2);
        assert_eq!(*p2.a.get(&x).unwrap(), 2.);
        assert_eq!(*p2.a.get(&y).unwrap(), 1.);
        assert!(!p2.a.contains_key(&z));

        let z3 = &x + &y + 40. + &z.cos();
        let p3 = z3.std_properties();
        assert!(!p3.affine);
        assert_eq!(p3.a.len(), 3);
        assert!(p3.a.contains_key(&x));
        assert!(p3.a.contains_key(&y));
        assert!(p3.a.contains_key(&z));
    }
}

