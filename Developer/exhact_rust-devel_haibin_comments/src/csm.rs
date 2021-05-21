/*!

Compressed sparse matrices (similar to CSC/CSR format)

# A sparse matrix oracle for compressed sparse row/column format

```
 Example 1: (to fill in later)
     - create a row-major csm with columns indexed by floats
     - access a row of the csm
     - access a column of the csm
```

*/




use std::iter::FromIterator;
use crate::matrix::{SmOracle, RingMetadata, MajorDimension};
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;

/**
A fairly conventional implementation of compressed sparse format (CSC and CSR) matrices, with a few differences:



* **Major dimension** is the dimension that's easiest to access.  
    * The **major fields** of a row-major `CSM` are its rows.  Row-major `CSM`s are almost the same as CSR's.
    * The **major fields** of a column-major `CSM` are its columns.  Column-major `CSM`s are almost the same as CSC's.

* **Minor keys** The major fields of a `CSM` are always indexed by integers. But the minor fields can be indexed by anything (e.g. simplices).  The `MinKey` type parameter records the type of the minor keys.

* **Coefficient ring** In order to accomodate a wide range of coefficient rings, each `CSM` stores some extra data about what its coefficient ring is and how to work with it.

The attributes of a `CSM` struct are as follows:

* `ringmetadata` = a struct that encodes all the information needed to work with the coefficient ring
* `nummaj` =  number of major indices (e.g., the number of rows if the `CSM` is row-major)
* `majdim` = major dimension (row or column)
* `snzval[majptr[i]..majptr[i+1]]` = the structurally nonzero coefficients of the `i`th major field 
* `snzind[majptr[i]..majptr[i+1]]` = the indices of the structurally nonzero coefficients of the `i`th major field 


```
Insert example here
```
*/
pub struct CSM<MinKey, SnzVal: Clone> {
    pub ringmetadata: RingMetadata<SnzVal>,
    pub nummaj: usize,                       // number of major indices
    pub majdim: MajorDimension,              // major dimension
    pub majptr: Vec<usize>,                  // maj end pointer
    pub minind: Vec<MinKey>,                 // min indices
    pub snzval: Vec<SnzVal>,                 // structural non-zero values
}

// Methods for CSM struct
impl<MinKey, SnzVal> CSM<MinKey, SnzVal> where
MinKey: Clone + Debug + Eq + PartialEq + Hash,
SnzVal: Clone + PartialEq + Debug
{
	/// Create a trivil CSM
    ///
    /// # Parameters
    /// - `majdim`: major dimension
    /// - `ringmetadata`: Ring meta data of structural non zero values in CSM
	pub fn new(
        majdim:             MajorDimension,
        ringmetadata:       RingMetadata<SnzVal>
    ) -> CSM<MinKey, SnzVal> {
		CSM {
            ringmetadata,
			nummaj: 0,
            majdim,
			majptr: vec![0],
			minind: Vec::new(),
			snzval: Vec::new()
		}
	}

    /// Create a trivil CSM with the specified capacity
    ///
    /// # Parameters
    /// - `capacity`: The capacity of allocation
    /// - `majdim`: major dimension
    /// - `ringmetadata`: Ring meta data of structural non zero values in CSM
    pub fn with_capacity(
        capacity: usize,
        majdim:             MajorDimension,
        ringmetadata:       RingMetadata<SnzVal>
    ) -> CSM<MinKey, SnzVal> {
        let mut majptr = Vec::with_capacity(capacity);
        majptr.push(0);
        CSM {
            ringmetadata,
            nummaj: 0,
            majdim,
            majptr,
            minind: Vec::with_capacity(capacity),
            snzval: Vec::with_capacity(capacity)
        }
    }

    /// Shrink the capacity of the CSM to save memory
    pub fn shrink_to_fit(&mut self) {
        self.majptr.shrink_to_fit();
        self.minind.shrink_to_fit();
        self.snzval.shrink_to_fit();
    }

    /// Print some fileds of CSM
    pub	fn print( &self ){
        println!("    Matrix in compressed format:");
        println!("    nummaj={};\n    majptr={:?}\n    minind={:?}\n    snzval={:?}.", self.nummaj, self.majptr, self.minind, self.snzval);
    }

	/// Reverse the order of majs
	pub fn reverse_maj_order(&mut self) {
		let mut reversed = CSM::new(self.majdim.clone(), self.ringmetadata.clone());
		for ii in 0..self.nummaj {
			let maj_index = self.nummaj-1-ii;
            reversed.append_maj(&mut self.maj_hash(&maj_index));
		}
		*self = reversed;
	}

	/// Add a new maj to the CSM
	///
	/// # Parameters
	/// - `hash`: A hash map representing the sparse major row
	pub fn append_maj(&mut self, hash: &mut HashMap<MinKey, SnzVal>) {
        for (key, val) in hash.drain() {
            self.push_snzval(key, val);
        }
		self.majptr.push(self.minind.len());
        self.nummaj += 1;
	}

	/// Add a single new entry to the CSM without updating the majptr
	///
	/// # Parameters
	/// - `ind`: The colum index of the entry
	/// - `val`: The value of the entry
	pub fn push_snzval(&mut self, ind:MinKey, val:SnzVal) { 
		self.minind.push(ind);
		self.snzval.push(val);
	}

    /// Change major dimension of a CSM matrix
    ///
    /// # Parameters
    /// - `minkey_list`: A list of minkeys that we pick as the majors of the new matrix
    pub fn change_major_dim(&self, minkey_list: &Vec<MinKey>) -> CSM<usize, SnzVal> {
        let mut output = CSM::new(MajorDimension::Col, self.ringmetadata.clone());
        for minkey in minkey_list.iter() {
            output.append_maj(&mut self.min_hash(minkey));
        }
        if self.majdim == MajorDimension::Col {
            output.majdim = MajorDimension::Row;
        }
        return output;
    }

}

/// An iterator that iterates all major keys corresponding to a given minor key
struct MajItr<'a, MinKey, SnzVal: Clone>{
    csm: &'a CSM<MinKey, SnzVal>,   // the underlying CSM matrix
    minind: MinKey,                 // the given minor key
    next_maj_ind: usize             // the next major key
}

/// implemenation of standard methods for MajItr struct
impl<MinInd: PartialEq, SnzVal: Clone> Iterator for MajItr<'_, MinInd, SnzVal>{
    type Item = (usize, SnzVal);

    /// Standard method of iterator to return the next item
    fn next( &mut self ) -> Option<Self::Item> {
        // look for the next row with a structurally nonzero value in this column
        for thisrow in self.next_maj_ind..self.csm.nummaj {
            for pointer in self.csm.majptr[thisrow]..self.csm.majptr[thisrow+1]{
                if self.csm.minind[pointer] == self.minind {
                    self.next_maj_ind = thisrow + 1;
                    return Some((thisrow, self.csm.snzval[pointer].clone()));
                }
            }
        }
        // if there is no such row, then return None and reset the iterator
        self.next_maj_ind = 0;
        return None;
    }

}

// See matrix.rs file for specific definition of SmOracle trait
impl< MinKey, SnzVal> SmOracle< usize, MinKey, SnzVal > for CSM< MinKey, SnzVal > where
MinKey: Clone + PartialEq + Eq + Hash,
SnzVal: Clone
{
    
	fn ring( &self ) -> &RingMetadata<SnzVal> {
        &self.ringmetadata
    }
	fn maj_dim( &self ) -> MajorDimension { MajorDimension::Row }

    fn maj_indmin( &self, majkey: &usize ) -> Option<Vec<MinKey>> {
        if self.majptr[*majkey]<self.majptr[*majkey+1] {
            return Some(Vec::from_iter(self.minind[self.majptr[*majkey]..self.majptr[*majkey+1]].iter().cloned()));
        } else {
            return None;
        }
    }

	fn maj_snzval( &self, majkey: &usize ) -> Option<Vec<SnzVal>> {
		if self.majptr[*majkey]<self.majptr[*majkey+1] {
			return Some(Vec::from_iter(self.snzval[self.majptr[*majkey]..self.majptr[*majkey+1]].iter().cloned()));
		} else {
			return None;
		}
	}

    fn min_indmaj( &self, minkey: &MinKey ) -> Option<Vec<usize>> {
        let mut colum = Vec::new();
        for ii in 0..(self.majptr.len()-1){
            for jj in self.majptr[ii]..self.majptr[ii+1]{
                if self.minind[jj] == *minkey {
                    colum.push(ii);
                }
            }
        }
        if colum.len() == 0 { return None; }
        else { return Some(colum); }
    }

    fn min_snzval( &self, minkey: &MinKey ) -> Option<Vec<SnzVal>> {
        let mut colum = Vec::new();
        for ii in 0..(self.majptr.len()-1){
            for jj in self.majptr[ii]..self.majptr[ii+1]{
                if self.minind[jj] == *minkey {
                    colum.push(self.snzval[jj].clone());
                }
            }
        }
        if colum.len() == 0 { return None; }
        else { return Some(colum); }
    }

	fn maj_itr( &self, majkey: &usize ) -> Box<dyn Iterator<Item=(MinKey, SnzVal)> + '_> {
		Box::new((self.majptr[*majkey]..self.majptr[*majkey+1]).map( move |x| (self.minind[x].clone(), self.snzval[x].clone())))
	}

	fn min_itr( &self, minkey: &MinKey ) -> Box<dyn Iterator<Item=(usize, SnzVal)> + '_> {
		Box::new( MajItr {
			csm: &self,
		    minind: minkey.clone(),
		    next_maj_ind: 0
		})
	}

    fn maj_fn( &self, majkey: &usize ) -> Box<dyn Fn(MinKey) -> Option<SnzVal> + '_> {
        let maj = *majkey;
		Box::new( move |x| {
				for ii in self.majptr[maj]..self.majptr[maj+1]{
					if self.minind[ii] == x {
						return Some(self.snzval[ii].clone());
					}
				}
				return None;
			}
		)
	}

	fn min_fn( &self, minkey: &MinKey ) -> Box<dyn Fn(usize) -> Option<SnzVal> + '_> {
        let min = minkey.clone();
		Box::new( move |x| {
				for ii in self.majptr[x]..self.majptr[x+1] {
					if self.minind[ii] == min {
						return Some(self.snzval[ii].clone());
					}
				}
				return None;
			}
		)
	}

	fn countsnz( &self ) -> Option<usize> {
		Some(self.snzval.len())
	}

	fn finiteminors( &self ) -> Option<bool> {
		Some(true)
	}
	fn finitemajors( &self ) -> Option<bool> {
		Some(true)
	}
}
