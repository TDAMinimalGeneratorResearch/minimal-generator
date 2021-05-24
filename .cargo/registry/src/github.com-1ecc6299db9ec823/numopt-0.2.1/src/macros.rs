//! Helper macros.

pub use approx;

/// Approximate equality of vectors based on
/// element-wise absolute difference.
#[macro_export]
macro_rules! assert_vec_approx_eq {
    ($x:expr, $y:expr, epsilon = $eps:expr) => {
        assert_eq!($x.len(), $y.len());
        for (a,b) in $x.iter().zip($y.iter()) {
            approx::assert_abs_diff_eq!(a, b, epsilon = $eps);
        }
    };
}
