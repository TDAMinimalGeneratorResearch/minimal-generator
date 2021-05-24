//! Trait for differentiating expression nodes.

use std::iter::FromIterator;
use std::collections::{HashMap, HashSet, VecDeque};

use crate::model::node::Node;
use crate::model::node_base::NodeBase;
use crate::model::constant::ConstantScalar;

/// Trait for differentiating expression nodes.
pub trait NodeDiff {

    /// Obtains all simple paths between node and given variable nodes.
    /// This is used for then performing differentiation using the chain-rule.
    fn all_simple_paths(&self, vars: &[&Node]) -> HashMap<Node, Vec<Vec<Node>>>;
    
    /// Obtains the derivative of the expression node with respect to a given variable node.
    fn derivative(&self, var: &Node) -> Node;

    /// Obtains the derivatives of the expression node with respect to a given array
    /// of variable nodes.
    fn derivatives(&self, vars: &[&Node]) -> HashMap<Node, Node>;
}

impl NodeDiff for Node {

    fn all_simple_paths(&self, vars: &[&Node]) -> HashMap<Node, Vec<Vec<Node>>> {

        // Check inputs
        for v in vars {
            match v {
                Node::VariableScalar(_x) => (),
                _ => panic!("variable expected")
            }
        }

        // Vars set
        let varset: HashSet<&Node> = HashSet::from_iter(vars.iter().map(|x| x.clone()));

        // Workqueue
        let mut wq: VecDeque<Vec<Node>> = VecDeque::new();
        wq.push_front(vec![self.clone()]);

        // Paths
        let mut paths: HashMap<Node, Vec<Vec<Node>>> = HashMap::new();
        for v in &varset {
            paths.insert((*v).clone(), Vec::new());
        }

        // Process
        loop {

            // Pop path
            let path = match wq.pop_front() {
                Some(p) => p,
                None => break
            };
            let node = path.last().unwrap();

            // Add paths
            match node {
                Node::VariableScalar(_x) => {
                    for v in &varset {
                        if node == *v {
                            let new_path = path.iter().map(|x| x.clone()).collect();
                            paths.get_mut(node).unwrap().push(new_path); 
                        }
                    } 
                }
                _ => (),
            };

            // Process arguments
            for n in node.arguments() {
                let mut new_path: Vec<Node> = path.iter().map(|x| x.clone()).collect();
                new_path.push(n.clone());
                wq.push_front(new_path);
            }
        }

        // Return paths
        paths
    }

    fn derivative(&self, var: &Node) -> Node {
        let derivs = self.derivatives(&vec![var]);
        derivs.get(var).unwrap().clone()
    }

    fn derivatives(&self, vars: &[&Node]) -> HashMap<Node, Node> {

        // Check inputs
        for v in vars {
            match v {
                Node::VariableScalar(_x) => (),
                _ => panic!("variable expected")
            }
        }

        // Vars set
        let varset: HashSet<&Node> = HashSet::from_iter(vars.iter().map(|x| x.clone()));

        // Derivatives
        let paths = self.all_simple_paths(vars);
        let mut derivs: HashMap<Node, Node> = HashMap::new();
        for v in varset.iter() {
            let mut d = ConstantScalar::new(0.);
            for path in paths.get(v).unwrap() {
                let mut prod = ConstantScalar::new(1.);
                for pair in path.as_slice().windows(2) {
                    prod = prod*pair[0].partial(&pair[1]);
                }
                d = d + prod;
            }
            derivs.insert((**v).clone(), d);
        }
        derivs
    }

}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::model::variable::VariableScalar;

    #[test]
    fn node_all_simple_paths() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let z = VariableScalar::new_continuous("z");

        let p1 = x.all_simple_paths(&vec![&x, &y, &z]);
        assert_eq!(p1.get(&x).unwrap().len(), 1);
        assert_eq!(p1.get(&x).unwrap()[0].len(), 1);
        assert_eq!(p1.get(&y).unwrap().len(), 0);
        assert_eq!(p1.get(&z).unwrap().len(), 0);

        let p2 = (&x + 1.).all_simple_paths(&vec![&x, &y, &z]);
        assert_eq!(p2.get(&x).unwrap().len(), 1);
        assert_eq!(p2.get(&x).unwrap()[0].len(), 2);
        assert_eq!(p2.get(&y).unwrap().len(), 0);
        assert_eq!(p2.get(&z).unwrap().len(), 0);

        let p3 = (4. + 3.*(&z + &x)).all_simple_paths(&vec![&x, &y, &z]);
        assert_eq!(p3.get(&x).unwrap().len(), 1);
        assert_eq!(p3.get(&x).unwrap()[0].len(), 3);
        assert_eq!(p3.get(&y).unwrap().len(), 0);
        assert_eq!(p3.get(&z).unwrap().len(), 1);
        assert_eq!(p3.get(&z).unwrap()[0].len(), 3);

        let f4 = &x + 5.;
        let g4 = &f4*(&z + 3.);
        let p4 = (f4 + g4).all_simple_paths(&vec![&x, &y, &z]);
        assert_eq!(p4.get(&x).unwrap().len(), 2);
        assert_eq!(p4.get(&x).unwrap()[0].len() +
                   p4.get(&x).unwrap()[1].len(), 6);
        assert_eq!(p4.get(&y).unwrap().len(), 0);
        assert_eq!(p4.get(&z).unwrap().len(), 1);
        assert_eq!(p4.get(&z).unwrap()[0].len(), 4);
    }
}