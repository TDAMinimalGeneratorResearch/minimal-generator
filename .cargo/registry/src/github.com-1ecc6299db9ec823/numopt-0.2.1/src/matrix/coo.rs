//! Sparse matrix in coordinate format.

use std::ops::Mul;

use crate::matrix::item::MatItem;
use crate::matrix::csr::CsrMat;

/// Sparse matrix in coordinate format.
#[derive(Debug, Clone)]
pub struct CooMat<T> {
    shape: (usize, usize),
    row_inds: Vec<usize>,
    col_inds: Vec<usize>,
    data: Vec<T>,
}

/// Iterator for nonzero elements of [CooMat](struct.CooMat.html)
pub struct CooMatIter<'a, T> {
    k: usize,
    mat: &'a CooMat<T>,
}

impl<T: MatItem> CooMat<T> {

    /// Creates [CooMat](struct.CooMat.html) from (row, col, value) triplets.
    pub fn new(shape: (usize, usize), 
               row_inds: Vec<usize>,
               col_inds: Vec<usize>,
               data: Vec<T>) -> Self {
        assert_eq!(row_inds.len(), col_inds.len());
        assert_eq!(row_inds.len(), data.len());
        Self {
            shape: shape,
            row_inds: row_inds,
            col_inds: col_inds,
            data: data,
        }
    }

    /// Creates [CooMat](struct.CooMat.html) from sparsity pattern.
    pub fn from_pattern(shape: (usize, usize), 
                        row_inds: Vec<usize>,
                        col_inds: Vec<usize>) -> Self {
        assert_eq!(row_inds.len(), col_inds.len());
        let data = vec![T::zero();row_inds.len()];
        Self {
            shape: shape,
            row_inds: row_inds,
            col_inds: col_inds,
            data: data,
        }
    }

    /// Creates empty [CooMat](struct.CooMat.html) from number of nonzero elements.
    pub fn from_nnz(shape: (usize, usize), nnz: usize) -> Self {
        Self {
            shape: shape,
            row_inds: vec![0;nnz],
            col_inds: vec![0;nnz],
            data: vec![T::zero();nnz],
        }
    }

    /// Number of rows.
    pub fn rows(&self) -> usize { self.shape.0 }

    /// Number of columns.
    pub fn cols(&self) -> usize { self.shape.1 }

    /// Number of nonzero elements.
    pub fn nnz(&self) -> usize { self.row_inds.len() }

    /// Vector of row indices.
    pub fn row_inds(&self) -> &[usize] { &self.row_inds }

    /// Vector of column indices.
    pub fn col_inds(&self) -> &[usize] { &self.col_inds }

    /// Vector of data values.
    pub fn data(&self) -> &[T] { &self.data }

    // Vector of mutable data values
    pub fn data_mut(&mut self) -> &mut [T] { &mut self.data }

    /// Sets row index.
    pub fn set_row_ind(&mut self, k: usize, row: usize) -> () { self.row_inds[k] = row }

    /// Sets column index.
    pub fn set_col_ind(&mut self, k: usize, row: usize) -> () { self.col_inds[k] = row }

    /// Creates iterator for nonzero elements.
    pub fn iter(&self) -> CooMatIter<T> { CooMatIter::new(&self) }

    /// Converts matrix to [CsrMat]((struct.CsrMat.html))
    pub fn to_csr(&self) -> CsrMat<T> {

        let mut indptr: Vec<usize> = vec![0; self.rows()+1];
        let mut indices: Vec<usize> = vec![0; self.nnz()];
        let mut data: Vec<T> = vec![T::zero(); self.nnz()];

        let mut counter: Vec<usize> = vec![0; self.rows()];

        // Count elements per row
        for row in self.row_inds.iter() {
            counter[*row] += 1;
        }
        
        // Set indptr
        indptr[0] = 0;
        let mut offset: usize = 0;
        for (i, c) in counter.iter().enumerate() {
            indptr[i+1] = offset + c;
            offset += c;
        }
        assert_eq!(indptr[self.rows()], self.nnz());

        // Set indices and data
        let mut k: usize; 
        counter.copy_from_slice(&vec![0; self.rows()]);
        for (row, col, val) in self.iter() {
            k = indptr[*row] + counter[*row]; 
            indices[k] = *col;
            data[k] = *val;
            counter[*row] += 1;
        }
        
        // Return
        CsrMat::new(
            self.shape,
            indptr,
            indices,
            data
        )
    }
}

impl<T: MatItem> Mul<Vec<T>> for &CooMat<T> {

    type Output = Vec<T>;

    fn mul(self, rhs: Vec<T>) -> Vec<T> {
        assert_eq!(self.cols(), rhs.len());
        let mut y = vec![T::zero(); self.rows()];
        for (row, col, val) in self.iter() {
            y[*row] += rhs[*col]*(*val);
        }
        y
    }
}

impl<'a, T: MatItem> CooMatIter<'a, T> {
    fn new(mat: &'a CooMat<T>) -> Self {
        Self {
            k: 0,
            mat: mat,
        }
    }
}

impl<'a, T: MatItem> Iterator for CooMatIter<'a, T> {
    type Item = (&'a usize, &'a usize, &'a T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.k < self.mat.nnz() {
            let item = (&self.mat.row_inds[self.k],
                        &self.mat.col_inds[self.k],
                        &self.mat.data[self.k]);
            self.k += 1;
            return Some(item);
        }
        else {
            self.k = 0;
            return None;
        }
    }
}

#[cfg(test)]
mod tests {
    
    use crate::matrix::coo::CooMat;
    use crate::assert_vec_approx_eq;

    #[test]
    fn coo_to_csr() {

        // 6 2 1 0 0
        // 3 1 0 7 0
        // 4 6 0 0 1

        let a = CooMat::new(
            (3, 5),
            vec![0 ,2 ,0 ,0 ,1 ,2  ,1 ,1 ,2 ,0 ,2],
            vec![0 ,1 ,2 ,0 ,0 ,4  ,1 ,3 ,0 ,1 ,4],
            vec![5.,6.,1.,1.,3.,-2.,1.,7.,4.,2.,3.],
        );

        let b = a.to_csr();

        assert_eq!(b.rows(), 3);
        assert_eq!(b.cols(), 5);
        assert_eq!(b.nnz(), 11);
        assert_vec_approx_eq!(b.indptr(),
                              vec![0, 4, 7, 11],
                              epsilon=0);
        assert_vec_approx_eq!(b.indices(),
                              vec![0, 2, 0, 1, 0, 1, 3, 1, 4, 0, 4],
                              epsilon=0);
        assert_vec_approx_eq!(b.data(),
                              vec![5., 1., 1., 2., 3., 1., 7., 6., -2., 4., 3.],
                              epsilon=1e-8);
    }

    #[test]
    fn coo_times_vec() {

        // 6 2 1 0 0
        // 3 1 0 7 0
        // 4 6 0 0 1

        let a = CooMat::new(
            (3, 5),
            vec![0 ,2 ,0 ,0 ,1 ,2  ,1 ,1 ,2 ,0 ,2],
            vec![0 ,1 ,2 ,0 ,0 ,4  ,1 ,3 ,0 ,1 ,4],
            vec![5.,6.,1.,1.,3.,-2.,1.,7.,4.,2.,3.],
        );

        let x = vec![2.,4.,3.,1.,7.];
        let y = (&a)*x;

        assert_vec_approx_eq!(y, vec![23., 17., 39.], epsilon = 1e-8);
    }
}