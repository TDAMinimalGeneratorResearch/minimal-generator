//! Optimization constraint.
//! 
//! Constraints currently supported are equality and inequality constraints
//! involving scalar expressions.

use std::fmt;
use std::ptr;
use std::rc::Rc;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use super::node::Node;
use super::node_base::NodeBase;

/// Constraint kind.
#[derive(PartialEq)]
pub enum ConstraintKind {
    Equal,
    LessEqual,
    GreaterEqual,
}

/// Constraint inner data.
struct ConstraintInner {
    lhs: Node,
    kind: ConstraintKind,
    rhs: Node,
    label: String,
}

/// Constraint.
pub struct Constraint(Rc<ConstraintInner>);

impl Constraint {

    /// Gets constraint kind.
    pub fn kind(&self) -> &ConstraintKind { &self.0.kind }

    /// Gets constraint label.
    pub fn label(&self) -> &str { self.0.label.as_ref() }

    /// Gets constraint left-hand-side.
    pub fn lhs(&self) -> &Node { &self.0.lhs }

    /// Creates new constraint.
    pub fn new(lhs: Node, kind: ConstraintKind, rhs: Node, label: &str) -> Constraint {
        Constraint(Rc::new(
            ConstraintInner{
                lhs: lhs,
                kind: kind,
                rhs: rhs,
                label: String::from(label),
            }
        ))
    }

    /// Gets constraint right-hand-side.
    pub fn rhs(&self) -> &Node { &self.0.rhs }

    /// Computes constraint violation given variable values.
    pub fn violation(&self, var_values: &HashMap<&Node, f64>) -> f64 {
        match self.0.kind {
            ConstraintKind::Equal => { 
                (self.0.lhs.evaluate(var_values)-self.0.rhs.evaluate(var_values)).abs()
            },
            ConstraintKind::LessEqual => { 
                0_f64.max(self.0.lhs.evaluate(var_values)-self.0.rhs.evaluate(var_values))
            },
            ConstraintKind::GreaterEqual => {
                0_f64.max(self.0.rhs.evaluate(var_values)-self.0.lhs.evaluate(var_values))
            },
        }
    }
}

impl Clone for Constraint {
    fn clone(&self) -> Self {
        Constraint(Rc::clone(&self.0))
    }
}

impl PartialEq for Constraint {

    fn eq(&self, other: &Self) -> bool {
        Rc::ptr_eq(&self.0, &other.0)
    }
}

impl Eq for ConstraintKind {}
impl Eq for Constraint {}

impl Hash for Constraint {
    
    fn hash<H: Hasher>(&self, state: &mut H) {
        ptr::hash(&*self.0, state);
    }
}

impl<'a> fmt::Display for Constraint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0.kind {
            ConstraintKind::Equal => write!(f, "{} == {}", self.0.lhs, self.0.rhs),
            ConstraintKind::LessEqual => write!(f, "{} <= {}", self.0.lhs, self.0.rhs),
            ConstraintKind::GreaterEqual => write!(f, "{} >= {}", self.0.lhs, self.0.rhs),
        }
    }
}

impl fmt::Debug for Constraint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self)
    }
}

#[cfg(test)]
mod tests {

    use maplit::hashmap;

    use super::*;
    use crate::model::variable::VariableScalar;
    use crate::model::constant::ConstantScalar;

    #[test]
    fn constr_hash() {

        let x = VariableScalar::new_continuous("x");

        let c = ConstantScalar::new(4.);

        let c1 = Constraint::new(x.clone(), ConstraintKind::Equal, c.clone(), "foo");
        let c2 = Constraint::new(x.clone(), ConstraintKind::LessEqual, c.clone(), "foo");
        let c3 = c1.clone();

        let h = hashmap!{ &c1 => 5., &c2 => 10.};
        assert_eq!(*h.get(&c1).unwrap(), 5.);
        assert_eq!(*h.get(&c2).unwrap(), 10.);
        assert_eq!(*h.get(&c3).unwrap(), 5.);
    }

    #[test]
    fn constr_clone_eq() {

        let x = VariableScalar::new_continuous("x");
        let c = ConstantScalar::new(4.);

        let c1 = Constraint::new(x.clone(), ConstraintKind::Equal, c.clone(), "foo");
        let c2 = Constraint::new(x.clone(), ConstraintKind::LessEqual, c.clone(), "foo");
        let c3 = c1.clone();

        assert_ne!(c1, c2);
        assert_eq!(c1, c3);
    }

    #[test]
    fn constr_label() {

        let x = VariableScalar::new_continuous("x");
        let c = ConstantScalar::new(4.);

        let z = Constraint::new(x, ConstraintKind::Equal, c, "foo");
        assert_eq!(z.label(), "foo");
    }

    #[test]
    fn constr_violation() {

        let x = VariableScalar::new_continuous("x");
        let c4 = ConstantScalar::new(4.);

        let var_values = hashmap!{ &x => 3. };

        let z1 = Constraint::new(x.clone(), ConstraintKind::Equal, c4.clone(), "foo");
        assert_eq!(z1.violation(&var_values), 1.);

        let z2 = Constraint::new(x.clone(), ConstraintKind::LessEqual, c4.clone(), "foo");
        assert_eq!(z2.violation(&var_values), 0.);

        let z3 = Constraint::new(x.clone(), ConstraintKind::LessEqual, -c4.clone(), "foo");
        assert_eq!(z3.violation(&var_values), 7.);

        let z4 = Constraint::new(x.clone(), ConstraintKind::GreaterEqual, c4.clone(), "foo");
        assert_eq!(z4.violation(&var_values), 1.);

        let z5 = Constraint::new(x.clone(), ConstraintKind::GreaterEqual, -c4.clone(), "foo");
        assert_eq!(z5.violation(&var_values), 0.);
    }
}


