//! Optimization algebraic modeling tools.
//!
//! These tools facilitate the construction, solution, and extraction of results
//! of optimization models.
//! 
//! ## Example
//! ```
//! use numopt::model::*;
//! use numopt::solver::*;
//! 
//! // Variables
//! let x = VariableScalar::new_continuous("x");
//! let y = VariableScalar::new_continuous("y");
//!
//! // Objective function
//! let f = 2.*&x - 5.*&y;
//! 
//! // Constraints
//! let c1 = x.geq(100.);
//! let c2 = x.leq(200.);
//! let c3 = y.geq(80.);
//! let c4 = y.leq(170.);
//! let c5 = (&y + 2.*&x).equal(&x + 200.);
//! 
//! // Model
//! let mut m = Model::new();
//! m.set_objective(Objective::minimize(&f));
//! m.add_constraint(&c1);
//! m.add_constraint(&c2);
//! m.add_constraint(&c3);
//! m.add_constraint(&c4);
//! m.add_constraint(&c5);
//! 
//! // Solver
//! let mut s = SolverClpCmd::new();
//! s.set_param("logLevel", SolverParam::IntParam(0)).unwrap();
//! m.solve(&s).unwrap();
//! 
//! // Status
//! assert_eq!(*m.solver_status().unwrap(), SolverStatus::Solved);
//! 
//! // Primal results
//! let final_primals = m.final_primals();
//! assert_eq!(*final_primals.get(&x).unwrap(), 100.);
//! assert_eq!(f.evaluate(&final_primals), -300.);
//! 
//! // Dual results
//! let final_duals = m.final_duals();
//! assert_eq!(*final_duals.get(&c1).unwrap(), 7.);
//! assert_eq!(*final_duals.get(&c5).unwrap(), -5.);
//! ```

pub mod node;
pub mod node_base;
pub mod node_func;
pub mod node_diff;
pub mod node_cmp;
pub mod node_std;
pub mod constant;
pub mod variable;
pub mod function;
pub mod constraint;
pub mod constraint_std;
pub mod model;
pub mod model_std;

pub use node::Node;
pub use node_cmp::NodeCmp;
pub use node_base::NodeBase;
pub use node_func::NodeFunc;
pub use node_diff::NodeDiff;
pub use variable::VariableScalar;
pub use constant::ConstantScalar;
pub use constraint::Constraint;
pub use model::Model;
pub use model::Objective;