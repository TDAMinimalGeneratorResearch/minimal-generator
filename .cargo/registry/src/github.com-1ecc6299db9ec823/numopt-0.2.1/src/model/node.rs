//! Optimization expression node.

use std::fmt;
use std::ptr;
use std::rc::Rc;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use num_traits::cast::ToPrimitive;
use num_traits::identities::{Zero, One};
use std::ops::{Add, Mul, Neg, Sub, Div};

use crate::model::node_base::NodeBase;
use crate::model::constant::ConstantScalar;
use crate::model::variable::VariableScalar;
use crate::model::function::add::FunctionAdd;
use crate::model::function::mul::FunctionMul;
use crate::model::function::div::FunctionDiv;
use crate::model::function::cos::FunctionCos;
use crate::model::function::sin::FunctionSin;

/// Expression node.
pub enum Node {
    ConstantScalar(Rc<ConstantScalar>),
    VariableScalar(Rc<VariableScalar>),
    FunctionAdd(Rc<FunctionAdd>),
    FunctionCos(Rc<FunctionCos>),
    FunctionDiv(Rc<FunctionDiv>),
    FunctionMul(Rc<FunctionMul>),
    FunctionSin(Rc<FunctionSin>),
}

impl Node {

    /// Determines whether node is a constant.
    pub fn is_constant(&self) -> bool {
        match self {
            Node::ConstantScalar(_) => true,
            _ => false
        }
    }

    /// Determines whether node is a constant with a given value.
    pub fn is_constant_with_value(&self, val: f64) -> bool {
        match self {
            Node::ConstantScalar(x) => x.value() == val,
            _ => false
        }
    }

    /// Gets node name. Currently, only variables have nonempty names.
    pub fn name(&self) -> &str {
        match self {
            Node::VariableScalar(x) => x.name(),
            _ => "",
        }
    }
}

impl Hash for Node {
    
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            Node::ConstantScalar(x) => ptr::hash(&**x, state),
            Node::VariableScalar(x) => ptr::hash(&**x, state),
            Node::FunctionAdd(x) => ptr::hash(&**x, state),
            Node::FunctionCos(x) => ptr::hash(&**x, state),
            Node::FunctionDiv(x) => ptr::hash(&**x, state),
            Node::FunctionMul(x) => ptr::hash(&**x, state),
            Node::FunctionSin(x) => ptr::hash(&**x, state),
        };
    }
}

impl PartialEq for Node {

    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Node::ConstantScalar(x), Node::ConstantScalar(y)) => Rc::ptr_eq(x, y),
            (Node::VariableScalar(x), Node::VariableScalar(y)) => Rc::ptr_eq(x, y),
            (Node::FunctionAdd(x), Node::FunctionAdd(y)) => Rc::ptr_eq(x, y),
            (Node::FunctionCos(x), Node::FunctionCos(y)) => Rc::ptr_eq(x, y),
            (Node::FunctionDiv(x), Node::FunctionDiv(y)) => Rc::ptr_eq(x, y),
            (Node::FunctionMul(x), Node::FunctionMul(y)) => Rc::ptr_eq(x, y),
            (Node::FunctionSin(x), Node::FunctionSin(y)) => Rc::ptr_eq(x, y),
            _ => false,
        }
    }
}

impl Eq for Node {}

impl Clone for Node {
    fn clone(&self) -> Self {
        match self {
            Node::ConstantScalar(x) => Node::ConstantScalar(Rc::clone(&x)),
            Node::VariableScalar(x) => Node::VariableScalar(Rc::clone(&x)),
            Node::FunctionAdd(x) => Node::FunctionAdd(Rc::clone(&x)),
            Node::FunctionCos(x) => Node::FunctionCos(Rc::clone(&x)),
            Node::FunctionDiv(x) => Node::FunctionDiv(Rc::clone(&x)),
            Node::FunctionMul(x) => Node::FunctionMul(Rc::clone(&x)),
            Node::FunctionSin(x) => Node::FunctionSin(Rc::clone(&x)), 
        }
    }
}

impl fmt::Display for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Node::ConstantScalar(x) => write!(f, "{}", x),
            Node::VariableScalar(x) => write!(f, "{}", x),
            Node::FunctionAdd(x) => write!(f, "{}", x),
            Node::FunctionCos(x) => write!(f, "{}", x),
            Node::FunctionDiv(x) => write!(f, "{}", x),
            Node::FunctionMul(x) => write!(f, "{}", x),
            Node::FunctionSin(x) => write!(f, "{}", x),
            
        }
    }
}

impl fmt::Debug for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Node::ConstantScalar(x) => write!(f, "{}", x),
            Node::VariableScalar(x) => write!(f, "{}", x),
            Node::FunctionAdd(x) => write!(f, "{}", x),
            Node::FunctionCos(x) => write!(f, "{}", x),
            Node::FunctionDiv(x) => write!(f, "{}", x),
            Node::FunctionMul(x) => write!(f, "{}", x),
            Node::FunctionSin(x) => write!(f, "{}", x),
        }
    }
}

macro_rules! impl_node_add_node {
    ($x: ty, $y: ty) => {
        impl Add<$y> for $x {
            type Output = Node;
            fn add(self, rhs: $y) -> Node {

                // Self zero
                if self.is_constant_with_value(0.) {
                    rhs.clone()
                }

                // Rhs zero
                else if rhs.is_constant_with_value(0.) {
                    self.clone()
                }

                // Other
                else {
                    let mut c: f64 = 0.;
                    let mut args: Vec<Node> = Vec::new();
                    for a in &[self.clone(), rhs.clone()] {
                        match a {
                            Node::ConstantScalar(x) => c += x.value(),
                            Node::FunctionAdd(x) => {
                                for x in x.arguments().iter() {
                                    match x {
                                        Node::ConstantScalar(y) => c += y.value(),
                                        _ => args.push((*x).clone()),
                                    }
                                }
                            },
                            _ => args.push(a.clone()),
                        };
                    }
                    if args.is_empty() {
                        return ConstantScalar::new(c);
                    }
                    if c != 0. {
                        args.push(ConstantScalar::new(c));
                    }
                    FunctionAdd::new(args)
                }
            }        
        }
    };
}

macro_rules! impl_node_add_scalar {
    ($x: ty, $y: ty) => {
        impl Add<$y> for $x {
            type Output = Node;
            fn add(self, rhs: $y) -> Node {

                // Rhs zero
                if rhs == <$y>::zero() {
                    self.clone()
                }

                // Other
                else {
                    let mut c: f64 = 0.;
                    let mut args: Vec<Node> = Vec::new();
                    for a in &[self.clone(), ConstantScalar::new(rhs.to_f64().unwrap())] {
                        match a {
                            Node::ConstantScalar(x) => c += x.value(),
                            Node::FunctionAdd(x) => {
                                for x in x.arguments().iter() {
                                    match x {
                                        Node::ConstantScalar(y) => c += y.value(),
                                        _ => args.push((*x).clone()),
                                    }
                                }
                            },
                            _ => args.push(a.clone()),
                        };        
                    }
                    if args.is_empty() {
                        return ConstantScalar::new(c);
                    }
                    if c != 0. {
                        args.push(ConstantScalar::new(c));
                    }
                    FunctionAdd::new(args)
                }
            }           
        }
        impl Add<$x> for $y {
            type Output = Node;
            fn add(self, rhs: $x) -> Node {

                // Self zero
                if self == <$y>::zero() {
                    rhs.clone()
                }

                // Other
                else {
                    let mut c: f64 = 0.;
                    let mut args: Vec<Node> = Vec::new();
                    for a in &[ConstantScalar::new(self.to_f64().unwrap()), rhs.clone()] {
                        match a {
                            Node::ConstantScalar(x) => c += x.value(),
                            Node::FunctionAdd(x) => {
                                for x in x.arguments().iter() {
                                    match x {
                                        Node::ConstantScalar(y) => c += y.value(),
                                        _ => args.push((*x).clone()),
                                    }
                                }
                            },
                            _ => args.push(a.clone()),
                        };    
                    }
                    if args.is_empty() {
                        return ConstantScalar::new(c);
                    }
                    if c != 0. {
                        args.push(ConstantScalar::new(c));
                    }
                    FunctionAdd::new(args)
                }
            }           
        }
    };
}

impl_node_add_node!(&Node, &Node);
impl_node_add_node!(&Node, Node);
impl_node_add_node!(Node, &Node);
impl_node_add_node!(Node, Node);
impl_node_add_scalar!(&Node, f64);
impl_node_add_scalar!(Node, f64);

macro_rules! impl_node_mul_node {
    ($x: ty, $y: ty) => {
        impl Mul<$y> for $x {
            type Output = Node;
            fn mul(self, rhs: $y) -> Node {

                // Self or rhs zero
                if self.is_constant_with_value(0.) || rhs.is_constant_with_value(0.) {
                    ConstantScalar::new(0.)
                }

                // Self one
                else if self.is_constant_with_value(1.) {
                    rhs.clone()
                }

                // Rhs one
                else if rhs.is_constant_with_value(1.) {
                    self.clone()
                }

                // Both constants
                else if self.is_constant() && rhs.is_constant() {
                    let vals: HashMap<&Node, f64> = HashMap::new();
                    ConstantScalar::new(self.evaluate(&vals)*rhs.evaluate(&vals))
                }

                // Other
                else {
                    let s = self.clone();
                    let r = rhs.clone();
                    match (&s, &r) {

                        // Constant times add
                        (Node::ConstantScalar(_x), Node::FunctionAdd(_y)) => {
                            FunctionAdd::new(r.arguments().iter().map(|x| &s*(*x)).collect())
                        },

                        // Add time constant
                        (Node::FunctionAdd(_x), Node::ConstantScalar(_y)) => {
                            FunctionAdd::new(s.arguments().iter().map(|x| (*x)*&r).collect())
                        },

                        // Other
                        _ => FunctionMul::new(s, r),
                    }
                }
            }        
        }
    };
}

macro_rules! impl_node_mul_scalar {
    ($x: ty, $y: ty) => {
        impl Mul<$y> for $x {
            type Output = Node;
            fn mul(self, rhs: $y) -> Node {

                // Self constant or rhs zero 
                if self.is_constant() || rhs == <$y>::zero() {
                    match self {
                        Node::ConstantScalar(x) => {
                            ConstantScalar::new(x.value()*rhs.to_f64().unwrap())
                        },
                        _ => ConstantScalar::new(0.),
                    }
                    
                }

                // Self one
                else if self.is_constant_with_value(1.) {
                    ConstantScalar::new(rhs.to_f64().unwrap())
                }

                // Rhs one
                else if rhs == <$y>::one() {
                    self.clone()
                }

                // Other
                else {
                    let s = self.clone();
                    let r = rhs.to_f64().unwrap();
                    match &s {

                        // Add times constant
                        Node::FunctionAdd(_x) => {
                            FunctionAdd::new(s.arguments().iter().map(|x| (*x)*r).collect())
                        },

                        // Other
                        _ => FunctionMul::new(s, ConstantScalar::new(r))
                    }
                }
            }           
        }
        impl Mul<$x> for $y {
            type Output = Node;
            fn mul(self, rhs: $x) -> Node {

                // Self zero or rhs constant 
                if self == <$y>::zero() || rhs.is_constant() {
                    match rhs {
                        Node::ConstantScalar(x) => {
                            ConstantScalar::new(self.to_f64().unwrap()*x.value())
                        },
                        _ => ConstantScalar::new(0.),
                    }
                }

                // Self one
                else if self == <$y>::one() {
                    rhs.clone()
                }

                // Rhs one
                else if rhs.is_constant_with_value(1.) {
                    ConstantScalar::new(self.to_f64().unwrap())
                }

                // Other
                else {
                    let s = self.to_f64().unwrap();
                    let r = rhs.clone();
                    match &r {

                        // Constant times add
                        Node::FunctionAdd(_x) => {
                            FunctionAdd::new(r.arguments().iter().map(|x| s*(*x)).collect())
                        },

                        // Other
                        _ => FunctionMul::new(ConstantScalar::new(s), r),
                    }
                }
            }           
        }
    };
}

impl_node_mul_node!(&Node, &Node);
impl_node_mul_node!(&Node, Node);
impl_node_mul_node!(Node, &Node);
impl_node_mul_node!(Node, Node);
impl_node_mul_scalar!(&Node, f64);
impl_node_mul_scalar!(Node, f64);

macro_rules! impl_node_neg {
    ($x: ty) => {
        impl Neg for $x {
            type Output = Node;
            fn neg(self) -> Node {
                match self {
                    Node::ConstantScalar(x) => {
                        ConstantScalar::new(-1.*x.value())
                    },
                    _ => (-1.)*self,
                }
            }        
        }
    };
}

impl_node_neg!(&Node);
impl_node_neg!(Node);

macro_rules! impl_node_sub_node {
    ($x: ty, $y: ty) => {
        impl Sub<$y> for $x {
            type Output = Node;
            fn sub(self, rhs: $y) -> Node {
                self + -1.*rhs
            }        
        }
    };
}

macro_rules! impl_node_sub_scalar {
    ($x: ty, $y: ty) => {
        impl Sub<$y> for $x {
            type Output = Node;
            fn sub(self, rhs: $y) -> Node {
                self + -1.*ConstantScalar::new(rhs.to_f64().unwrap())
            }           
        }
        impl Sub<$x> for $y {
            type Output = Node;
            fn sub(self, rhs: $x) -> Node {
                ConstantScalar::new(self.to_f64().unwrap()) + -1.*rhs
            }           
        }
    };
}

impl_node_sub_node!(&Node, &Node);
impl_node_sub_node!(&Node, Node);
impl_node_sub_node!(Node, &Node);
impl_node_sub_node!(Node, Node);
impl_node_sub_scalar!(&Node, f64);
impl_node_sub_scalar!(Node, f64);

macro_rules! impl_node_div_node {
    ($x: ty, $y: ty) => {
        impl Div<$y> for $x {
            type Output = Node;
            fn div(self, rhs: $y) -> Node {

                // Rhs is zero
                if rhs.is_constant_with_value(0.) {
                    panic!("dividion by zero constant")
                }

                // Rhs is one
                else if rhs.is_constant_with_value(1.) {
                    self.clone()
                }

                // Self is zero
                else if self.is_constant_with_value(0.) {
                    ConstantScalar::new(0.)
                }

                // Both are constants
                else if self.is_constant() && rhs.is_constant() {
                    let vals: HashMap<&Node, f64> = HashMap::new();
                    ConstantScalar::new(self.evaluate(&vals)/rhs.evaluate(&vals))
                }

                // Other
                else {
                    FunctionDiv::new(self.clone(), rhs.clone())
                }
            }        
        }
    };
}

macro_rules! impl_node_div_scalar {
    ($x: ty, $y: ty) => {
        impl Div<$y> for $x {
            type Output = Node;
            fn div(self, rhs: $y) -> Node {

                // Rhs is zero
                if rhs == <$y>::zero() {
                    panic!("dividion by zero constant")
                }

                // Rhs is one
                else if rhs == <$y>::one() {
                    self.clone()
                }

                // Self is zero
                else if self.is_constant_with_value(0.) {
                    ConstantScalar::new(0.)
                }

                // Both are constants
                else if self.is_constant() {
                    let vals: HashMap<&Node, f64> = HashMap::new();
                    ConstantScalar::new(self.evaluate(&vals)/rhs.to_f64().unwrap())
                }

                // Other
                else {
                    FunctionDiv::new(self.clone(), 
                                     ConstantScalar::new(rhs.to_f64().unwrap()))
                }
            }           
        }
        impl Div<$x> for $y {
            type Output = Node;
            fn div(self, rhs: $x) -> Node {

                // Rhs is zero
                if rhs.is_constant_with_value(0.) {
                    panic!("dividion by zero constant")
                }

                // Rhs is oen
                else if rhs.is_constant_with_value(1.) {
                    ConstantScalar::new(self.to_f64().unwrap())
                }

                // Self is zero
                else if self == <$y>::zero() {
                    ConstantScalar::new(0.)
                }

                // Both are constants
                else if rhs.is_constant() {
                    let vals: HashMap<&Node, f64> = HashMap::new();
                    ConstantScalar::new(self.to_f64().unwrap()/rhs.evaluate(&vals))
                }

                // Other
                else {
                    FunctionDiv::new(ConstantScalar::new(self.to_f64().unwrap()), 
                                     rhs.clone())
                }
            }           
        }
    };
}

impl_node_div_node!(&Node, &Node);
impl_node_div_node!(&Node, Node);
impl_node_div_node!(Node, &Node);
impl_node_div_node!(Node, Node);
impl_node_div_scalar!(&Node, f64);
impl_node_div_scalar!(Node, f64);

#[cfg(test)]
mod tests {

    use maplit::hashmap;
    use num_traits::pow::Pow;

    use crate::model::node_base::NodeBase;
    use crate::model::variable::VariableScalar;
    use crate::model::constant::ConstantScalar;

    #[test]
    fn node_hash() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let xx = x.clone();

        let h = hashmap!{ &x => 2., &y => 4.};
        assert_eq!(*h.get(&x).unwrap(), 2.);
        assert_eq!(*h.get(&y).unwrap(), 4.);
        assert_eq!(*h.get(&xx).unwrap(), 2.);
    }

    #[test]
    fn node_add_node() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let c0 = ConstantScalar::new(0.);
        let c1 = ConstantScalar::new(1.);
        let c2 = ConstantScalar::new(2.);

        let var_values = hashmap!{ &x => 3., &y => 4. };

        let z1 = &x + &y;
        assert_eq!(format!("{}", z1), "x + y");
        assert_eq!(z1.evaluate(&var_values), 7.);

        let z2 = &y + &x;
        assert_eq!(format!("{}", z2), "y + x");
        assert_eq!(z2.evaluate(&var_values), 7.);

        let z3 = &x + (&y + &x);
        assert_eq!(format!("{}", z3), "x + y + x");
        assert_eq!(z3.evaluate(&var_values), 10.);

        let z4 = (&x + &y) + &x;
        assert_eq!(format!("{}", z4), "x + y + x");
        assert_eq!(z4.evaluate(&var_values), 10.);

        let z5 = &z1 + &z2 + &z3 + &z4;
        assert_eq!(format!("{}", z5), "x + y + y + x + x + y + x + x + y + x");
        assert_eq!(z5.evaluate(&var_values), 34.);
    
        let z6 = (&x + &y) + (&y + &x);
        assert_eq!(format!("{}", z6), "x + y + y + x");
        assert_eq!(z6.evaluate(&var_values), 14.);

        let z7 = &x + &c0;
        assert_eq!(z7, x);

        let z8 = &c0 + &x;
        assert_eq!(z8, x);

        let z9 = (&x + 1.) + &y;
        assert_eq!(format!("{:?}", z9.arguments()), "[x, y, 1]");
        assert_eq!(z9.evaluate(&var_values), 8.);

        let z10 = &x + (&y + 5.);
        assert_eq!(format!("{:?}", z10.arguments()), "[x, y, 5]");
        assert_eq!(z10.evaluate(&var_values), 12.);

        let z11 = (&x + 2.) + (&y + 7.);
        assert_eq!(format!("{:?}", z11.arguments()), "[x, y, 9]");
        assert_eq!(z11.evaluate(&var_values), 16.);

        let z12 = &c1 + &c2;
        assert!(z12.is_constant_with_value(3.));

        let z13 = (&x + 4.) + (5. + &y);
        assert_eq!(format!("{}", z13), "x + y + 9");
        assert_eq!(z13.evaluate(&var_values), 16.);
    }

    #[test]
    fn node_add_scalar() {

        let x = VariableScalar::new_continuous("x");
        let c1 = ConstantScalar::new(1.);

        let var_values = hashmap!{ &x => 3. };

        let z1 = &x + 15.;
        assert_eq!(format!("{}", z1), "x + 15");
        assert_eq!(z1.evaluate(&var_values), 18.);

        let z2 = 13. + &x;
        assert_eq!(format!("{}", z2), "x + 13");
        assert_eq!(z2.evaluate(&var_values), 16.);

        let z3 = 2. + &z2 + 6.;
        assert_eq!(format!("{}", z3), "x + 21");
        assert_eq!(z3.evaluate(&var_values), 24.);

        let z4 = &x + 0.;
        assert_eq!(z4, x);

        let z5 = 0. + &x;
        assert_eq!(z5, x);

        let z6 = (&x + 1.) + 2.;
        assert_eq!(format!("{:?}", z6.arguments()), "[x, 3]");
        assert_eq!(z6.evaluate(&var_values), 6.);

        let z7 = 3. + (&x + 4.);
        assert_eq!(format!("{:?}", z7.arguments()), "[x, 7]");
        assert_eq!(z7.evaluate(&var_values), 10.);

        let z8 = 4. + &c1;
        assert!(z8.is_constant_with_value(5.));

        let z9 = &c1 + 5.;
        assert!(z9.is_constant_with_value(6.));

        let z10 = (&x + 4.) + 5.;
        assert_eq!(format!("{}", z10), "x + 9");
        assert_eq!(z10.evaluate(&var_values), 12.);
        
        let z11 = 3. + (&x + 4.) + 5.;
        assert_eq!(format!("{}", z11), "x + 12");
        assert_eq!(z11.evaluate(&var_values), 15.);
    }

    #[test]
    fn node_mul_node() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let c0 = ConstantScalar::new(0.);
        let c1 = ConstantScalar::new(1.);
        let c2 = ConstantScalar::new(2.);
        let c3 = ConstantScalar::new(3.);

        let var_values = hashmap!{ &x => 3., &y => 4. };

        let z1 = &x*&y;
        assert_eq!(format!("{}", z1), "x*y");
        assert_eq!(z1.evaluate(&var_values), 12.);

        let z2 = &y*&x;
        assert_eq!(format!("{}", z2), "y*x");
        assert_eq!(z2.evaluate(&var_values), 12.);

        let z3 = (&y*&x)*&x;
        assert_eq!(format!("{}", z3), "y*x*x");
        assert_eq!(z3.evaluate(&var_values), 36.);

        let z4 = &y*(&x*&x);
        assert_eq!(format!("{}", z4), "y*x*x");
        assert_eq!(z4.evaluate(&var_values), 36.);

        let z5 = &z4*(&x*&z3);
        assert_eq!(format!("{}", z5), "y*x*x*x*y*x*x");
        assert_eq!(z5.evaluate(&var_values), (4.).pow(2.)*((3.).pow(5.)));

        let z6 = &x*&c0;
        assert!(z6.is_constant_with_value(0.));

        let z7 = &c0*&x;
        assert!(z7.is_constant_with_value(0.));

        let z8 = &c1*&x;
        assert_eq!(z8, x);

        let z9 = &x*&c1;
        assert_eq!(z9, x);

        let z10 = &c1*&c0;
        assert!(z10.is_constant_with_value(0.));

        let z11 = &c0*&c1;
        assert!(z11.is_constant_with_value(0.));

        let z12 = (&x + &y*&x)*(&y*&x + &y);
        assert_eq!(format!("{}", z12), "(x + y*x)*(y*x + y)");
        assert_eq!(z12.evaluate(&var_values), 15.*16.);

        let z13 = &c3*&c2;
        assert!(z13.is_constant_with_value(6.));

        let z14 = &c3*(&x + 3.);
        assert_eq!(format!("{}", z14), "3*x + 9");
        assert_eq!(z14.evaluate(&var_values), 18.);

        let z15 = (&x + &y)*&c2;
        assert_eq!(format!("{}", z15), "x*2 + y*2");
        assert_eq!(z15.evaluate(&var_values), 14.);
    }

    #[test]
    fn node_mul_scalar() {

        let x = VariableScalar::new_continuous("x");
        let c1 = ConstantScalar::new(1.);
        let c2 = ConstantScalar::new(2.);

        let var_values = hashmap!{ &x => 3. };

        let z1 = &x*15.;
        assert_eq!(format!("{}", z1), "x*15");
        assert_eq!(z1.evaluate(&var_values), 45.);

        let z2 = 13.*&x;
        assert_eq!(format!("{}", z2), "13*x");
        assert_eq!(z2.evaluate(&var_values), 39.);

        let z3 = 2.*&z2*6.;
        assert_eq!(format!("{}", z3), "2*13*x*6");
        assert_eq!(z3.evaluate(&var_values), 2.*13.*3.*6.);

        let z4 = &x*0.;
        assert!(z4.is_constant_with_value(0.));

        let z5 = 0.*&x;
        assert!(z5.is_constant_with_value(0.));

        let z6 = &x*1.;
        assert_eq!(z6, x);

        let z7 = 1.*&x;
        assert_eq!(z7, x);

        let z8 = &c1*0.;
        assert!(z8.is_constant_with_value(0.));

        let z9 = 0.*&c1;
        assert!(z9.is_constant_with_value(0.));

        let z10 = 3.*&c2;
        assert!(z10.is_constant_with_value(6.));

        let z11 = &c2*5.;
        assert!(z11.is_constant_with_value(10.));

        let z12 = 4.*(&x + 3.);
        assert_eq!(format!("{}", z12), "4*x + 12");
        assert_eq!(z12.evaluate(&var_values), 24.);

        let z13 = (4. + &x)*10.;
        assert_eq!(format!("{}", z13), "x*10 + 40");
        assert_eq!(z13.evaluate(&var_values), 70.);
    }

    #[test]
    fn node_neg() {

        let x = VariableScalar::new_continuous("x");
        let c = ConstantScalar::new(5.);

        let var_values = hashmap!{ &x => 3. };

        let z1 = -&x;
        assert_eq!(format!("{}", z1), "-1*x");
        assert_eq!(z1.evaluate(&var_values), -3.);

        let z2 = -(&x + 3.);
        assert_eq!(format!("{}", z2), "-1*x + -3");
        assert_eq!(z2.evaluate(&var_values), -6.);

        let z3 = -&c;
        assert!(z3.is_constant_with_value(-5.));
    }

    #[test]
    fn node_sub_node() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");

        let var_values = hashmap!{ &x => 3., &y => 4. };

        let z1 = &x - &y;
        assert_eq!(z1.evaluate(&var_values), -1.);
        assert_eq!(format!("{}", z1), "x + -1*y");

        let z2 = &y - &x;
        assert_eq!(z2.evaluate(&var_values), 1.);
        assert_eq!(format!("{}", z2), "y + -1*x");

        let z3 = &x - (&x - &y);
        assert_eq!(z3.evaluate(&var_values), 4.);
        assert_eq!(format!("{}", z3), "x + -1*x + -1*-1*y");

        let z4 = (&x - &y) - &y;
        assert_eq!(z4.evaluate(&var_values), -5.);

        let z5 = &z4 - &z3 - &x;
        assert_eq!(z5.evaluate(&var_values), -12.);

        let z6 = (&z1 - &z2) - (&z3 - &z4);
        assert_eq!(z6.evaluate(&var_values), -2.-9.);
    }

    #[test]
    fn node_sub_scalar() {

        let x = VariableScalar::new_continuous("x");

        let var_values = hashmap!{ &x => 3. };

        let z1 = &x - 15.;
        assert_eq!(format!("{}", z1), "x + -15");
        assert_eq!(z1.evaluate(&var_values), -12.);

        let z2 = 13. - &x;
        assert_eq!(format!("{}", z2), "-1*x + 13");
        assert_eq!(z2.evaluate(&var_values), 10.);

        let z3 = 2. - &z2 - 6.;
        assert_eq!(format!("{}", z3), "-1*-1*x + -17");
        assert_eq!(z3.evaluate(&var_values), -14.);
    }

    #[test]
    fn node_div_node() {

        let x = VariableScalar::new_continuous("x");
        let y = VariableScalar::new_continuous("y");
        let c0 = ConstantScalar::new(0.);
        let c1 = ConstantScalar::new(1.);
        let c2 = ConstantScalar::new(2.);
        let c3 = ConstantScalar::new(3.);

        let var_values = hashmap!{ &x => 3., &y => 4. };

        let z1 = &x/&y;
        assert_eq!(format!("{}", z1), "x/y");
        assert_eq!(z1.evaluate(&var_values), 3./4.);

        let z2 = (3.*&x)/(4.*&y);
        assert_eq!(format!("{}", z2), "3*x/(4*y)");
        assert_eq!(z2.evaluate(&var_values), 9./16.);

        let z3 = (3. + &x)/(&y + 4.);
        assert_eq!(format!("{}", z3), "(x + 3)/(y + 4)");
        assert_eq!(z3.evaluate(&var_values), 6./8.);

        let z4 = &x/(3.+&y);
        assert_eq!(format!("{}", z4), "x/(y + 3)");
        assert_eq!(z4.evaluate(&var_values), 3./7.);

        let z5 = (2.+&x)/&y;
        assert_eq!(format!("{}", z5), "(x + 2)/y");
        assert_eq!(z5.evaluate(&var_values), 5./4.);

        let z6 = &x/&c1;
        assert_eq!(z6, x);

        let z7 = &c0/&x;
        assert!(z7.is_constant_with_value(0.));

        let z8 = &c2/&c3;
        assert!(z8.is_constant_with_value(2./3.));
    }

    #[test]
    fn node_div_scalar() {

        let x = VariableScalar::new_continuous("x");
        let c1 = ConstantScalar::new(1.);
        let c2 = ConstantScalar::new(2.);

        let var_values = hashmap!{ &x => 4. };

        let z1 = 3./&x;
        assert_eq!(format!("{}", z1), "3/x");
        assert_eq!(z1.evaluate(&var_values), 3./4.);

        let z2 = 3./(&x + 1.);
        assert_eq!(format!("{}", z2), "3/(x + 1)");
        assert_eq!(z2.evaluate(&var_values), 3./5.);

        let z3 = &x/3.;
        assert_eq!(format!("{}", z3), "x/3");
        assert_eq!(z3.evaluate(&var_values), 4./3.);

        let z4 = (&x + 1.)/3.;
        assert_eq!(format!("{}", z4), "(x + 1)/3");
        assert_eq!(z4.evaluate(&var_values), 5./3.);

        let z5 = &x/1.;
        assert_eq!(z5, x);

        let z6 = 0./&x;
        assert!(z6.is_constant_with_value(0.));

        let z7 = 0./&c1;
        assert!(z7.is_constant_with_value(0.));

        let z8 = &c2/4.;
        assert!(z8.is_constant_with_value(2./4.));

        let z9 = 5./&c2;
        assert!(z9.is_constant_with_value(5./2.));
    }
}