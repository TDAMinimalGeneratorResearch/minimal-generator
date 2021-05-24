use super::Itertools;
use super::size_hint;
use super::misc::GenericRange;

/// An iterator adaptor that yields a subset (a slice) of the base iterator.
///
/// **Note:** slicing an iterator is not constant time, and much less efficient than
/// slicing for example a vector.
///
/// ```
/// use std::iter::repeat;
/// use itertools::Itertools;
///
/// let it = repeat('a').slice(..3);
/// assert_eq!(it.count(), 3);
/// ```
///
/// See [`.slice()`](trait.Itertools.html#method.slice) for more information.
#[derive(Copy, Clone)]
pub struct ISlice<I> {
    start: usize,
    end: usize,
    iter: I,
}

impl<I> ISlice<I>
    where I: Iterator
{
    /// Create a new **ISlice**.
    pub fn new<R: GenericRange>(iter: I, range: R) -> Self {
        let mut start = range.start().unwrap_or(0);
        let end = range.end().unwrap_or(::std::usize::MAX);
        if start > end {
            start = end;
        }
        ISlice {
            start: start,
            end: end,
            iter: iter,
        }
    }
}

impl<I> Iterator for ISlice<I>
    where I: Iterator
{
    type Item = I::Item;

    fn next(&mut self) -> Option<I::Item> {
        if self.start != 0 {
            let st = self.start;
            let n = self.iter.dropn(self.start);
            self.start = 0;
            self.end -= n;
            if n != st {
                // iterator is already done.
                return None;
            }
        }
        if self.end != 0 {
            self.end -= 1;
            self.iter.next()
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.end - self.start;
        let sh = self.iter.size_hint();
        size_hint::min(size_hint::sub_scalar(sh, self.start), (len, Some(len)))
    }
}

impl<I> ExactSizeIterator for ISlice<I>
    where I: ExactSizeIterator
{}
