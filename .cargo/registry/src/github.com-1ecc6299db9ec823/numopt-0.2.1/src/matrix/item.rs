//! Sparse matrix element.

use num_traits::{Zero};
use std::ops::{Mul, AddAssign};

/// A trait for sparse matrix elements.
pub trait MatItem: Zero + Mul<Output=Self> + AddAssign + Clone + Copy {}

impl<T: Zero + Mul<Output=Self> + AddAssign + Clone + Copy> MatItem for T {}