//! Sparse matrix in compressed sparse row format.

use crate::matrix::item::MatItem;

/// Sparse matrix in compressed sparse row format.
#[derive(Debug, Clone)]
pub struct CsrMat<T> {
    shape: (usize, usize),
    indptr: Vec<usize>,
    indices: Vec<usize>,
    data: Vec<T>,
}

impl<T: MatItem> CsrMat<T> {

    /// Creates [CsrMat](struct.CsrMat.html) from raw data.
    pub fn new(shape: (usize, usize), 
               indptr: Vec<usize>,
               indices: Vec<usize>,
               data: Vec<T>) -> Self {
        assert_eq!(indptr.len(), shape.0+1);
        assert_eq!(indices.len(), data.len());
        assert_eq!(*indptr.last().unwrap(), data.len());
        Self {
            shape: shape,
            indptr: indptr,
            indices: indices,
            data: data,
        }
    }

    /// Number of rows.
    pub fn rows(&self) -> usize { self.shape.0 }

    /// Number of columns.
    pub fn cols(&self) -> usize { self.shape.1 }

    /// Number of nonzero elements.
    pub fn nnz(&self) -> usize { self.indices.len() }

    /// Vector of index pointers.
    pub fn indptr(&self) -> &[usize] { &self.indptr }

    /// Vector of column indices.
    pub fn indices(&self) -> &[usize] { &self.indices }

    /// Vector of data values.
    pub fn data(&self) -> &[T] { &self.data }

    /// Sums duplicate entries in-place.
    pub fn sum_duplicates(&mut self) -> () {

        let mut colseen: Vec<bool> = vec![false; self.cols()];
        let mut colrow: Vec<usize> = vec![0; self.cols()];
        let mut colnewk: Vec<usize> = vec![0; self.cols()];

        let mut d: T;
        let mut col: usize;
        let mut new_k: usize = 0;
        let mut new_counter: Vec<usize> = vec![0; self.rows()];
        let mut new_indices: Vec<usize> = Vec::new();
        let mut new_data: Vec<T> = Vec::new();
        for row in 0..self.rows() {
            for k in self.indptr[row]..self.indptr[row+1] {
                
                col = self.indices[k];
                d = self.data[k];

                // New column in row
                if !colseen[col] || colrow[col] != row {        
                    colnewk[col] = new_k;
                    new_counter[row] += 1;
                    new_indices.push(col);
                    new_data.push(d);
                    new_k += 1;
                }
                
                // Duplicate column in row
                else { 
                    new_data[colnewk[col]] += d;
                }

                // Update
                colseen[col] = true;
                colrow[col] = row;
            }

        }

        let mut offset: usize = 0;
        let mut new_indptr: Vec<usize> = vec![0; self.rows()+1];
        for (row, c) in new_counter.iter().enumerate() {
            new_indptr[row+1] = offset + c;
            offset += c;
        }

        self.indptr = new_indptr;
        self.indices = new_indices;
        self.data = new_data;

        assert_eq!(self.indptr.len(), self.rows()+1);
        assert_eq!(self.indices.len(), self.indptr[self.rows()]);
        assert_eq!(self.indices.len(), self.data.len());
    }
}

#[cfg(test)]
mod tests {

    use crate::matrix::coo::CooMat;
    use crate::assert_vec_approx_eq;

    #[test]
    fn csr_sum_dublicates() {

        // 6 2 1 0 0
        // 3 1 0 7 0
        // 4 6 0 0 1

        let a = CooMat::new(
            (3, 5),
            vec![0 ,2 ,0 ,0 ,1 ,2  ,1 ,1 ,2 ,0 ,2],
            vec![0 ,1 ,2 ,0 ,0 ,4  ,1 ,3 ,0 ,1 ,4],
            vec![5.,6.,1.,1.,3.,-2.,1.,7.,4.,2.,3.],
        );

        let mut b = a.to_csr();
        b.sum_duplicates();

        assert_eq!(b.rows(), 3);
        assert_eq!(b.cols(), 5);
        assert_eq!(b.nnz(), 9);
        assert_vec_approx_eq!(b.indptr(),
                              vec![0, 3, 6, 9],
                              epsilon=0);
        assert_vec_approx_eq!(b.indices(),
                              vec![0, 2, 1, 0, 1, 3, 1, 4, 0],
                              epsilon=0);
        assert_vec_approx_eq!(b.data(),
                              vec![6., 1., 2., 3., 1., 7., 6., 1., 4.],
                              epsilon=1e-8);
    }
}

