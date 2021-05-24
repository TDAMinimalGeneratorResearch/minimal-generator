//! This module allows solving problems is external solver binaries.
//! Contrarily to other solver modules, this one doesn't require linking your program to any solver.
//! A solver binary will need to be present on the user's computer at runtime.

use std::cmp::Ordering;

use lp_solvers::lp_format::LpObjective;
use lp_solvers::problem::StrExpression;
pub use lp_solvers::solvers::*;
use lp_solvers::util::UniqueNameGenerator;

use crate::constraint::ConstraintReference;
use crate::solvers::ObjectiveDirection;
use crate::variable::UnsolvedProblem;
use crate::{
    Constraint, Expression, IntoAffineExpression, ResolutionError, Solution, Solver, SolverModel,
    Variable,
};

/// An external solver
pub struct LpSolver<T: lp_solvers::solvers::SolverTrait>(pub T);

impl<T: lp_solvers::solvers::SolverTrait + Clone> Solver for LpSolver<T> {
    type Model = Model<T>;

    fn create_model(&mut self, problem: UnsolvedProblem) -> Self::Model {
        let name = "good_lp_problem".to_string();
        let sense = match problem.direction {
            ObjectiveDirection::Maximisation => LpObjective::Maximize,
            ObjectiveDirection::Minimisation => LpObjective::Minimize,
        };
        let mut gen = UniqueNameGenerator::default();
        let variables: Vec<lp_solvers::problem::Variable> = problem
            .variables
            .into_iter()
            .map(|def| lp_solvers::problem::Variable {
                name: gen.add_variable(&def.name).to_string(),
                is_integer: def.is_integer,
                lower_bound: def.min,
                upper_bound: def.max,
            })
            .collect();
        let objective = linear_coefficients_str(&problem.objective, &variables);
        Model {
            problem: lp_solvers::problem::Problem {
                name,
                sense,
                objective,
                variables,
                constraints: vec![],
            },
            solver: self.0.clone(),
        }
    }
}

/// A problem to be used by lp-solvers
pub struct Model<T> {
    problem: lp_solvers::problem::Problem,
    solver: T,
}

impl<T: SolverTrait> SolverModel for Model<T> {
    type Solution = LpSolution;
    type Error = ResolutionError;

    fn solve(self) -> Result<Self::Solution, Self::Error> {
        let map = self.solver.run(&self.problem)?;
        match map.status {
            Status::Infeasible => return Err(ResolutionError::Infeasible),
            Status::Unbounded => return Err(ResolutionError::Unbounded),
            Status::NotSolved => return Err(ResolutionError::Other("unknown error: not solved")),
            _ => {}
        }
        let solution = self
            .problem
            .variables
            .iter()
            .map(|v| f64::from(*map.results.get(&v.name).unwrap_or(&0.)))
            .collect();
        Ok(LpSolution { solution })
    }

    fn add_constraint(&mut self, c: Constraint) -> ConstraintReference {
        let reference = ConstraintReference {
            index: self.problem.constraints.len(),
        };
        self.problem
            .constraints
            .push(lp_solvers::lp_format::Constraint {
                lhs: linear_coefficients_str(&c.expression, &self.problem.variables),
                operator: if c.is_equality {
                    Ordering::Equal
                } else {
                    Ordering::Less
                },
                rhs: -c.expression.constant,
            });
        reference
    }
}

fn linear_coefficients_str(
    expr: &Expression,
    variables: &[lp_solvers::problem::Variable],
) -> StrExpression {
    StrExpression(
        expr.linear_coefficients()
            .map(|(var, coeff)| format!("{} {}", coeff, variables[var.index()].name))
            .collect::<Vec<String>>()
            .join(" + "),
    )
}

/// A solution
pub struct LpSolution {
    solution: Vec<f64>,
}

impl Solution for LpSolution {
    fn value(&self, variable: Variable) -> f64 {
        self.solution[variable.index()]
    }
}
