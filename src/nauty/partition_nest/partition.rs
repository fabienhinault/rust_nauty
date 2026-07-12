use crate::nauty::Set;

use super::PartitionNest;
use super::partition_nest_chunk_by::PartitionNestChunkBy;
use std::ops::Index;
use std::ops::IndexMut;

#[derive(Default, PartialEq, Debug)]
pub struct Partition {
    pub(crate) nest: PartitionNest,
    pub(crate) level: usize,
}

impl Partition {
    pub fn new(mut nest: PartitionNest, level: usize) -> Self {
        nest.extend_numcells(level);
        Self { nest, level }
    }

    pub fn numcells(&self) -> usize {
        self.nest.numcells(self.level)
    }

    pub fn numcells_mut(&mut self) -> &mut usize {
        self.nest.numcells_mut(self.level)
    }

    /// last index of cell containing i for partition of given level
    pub fn cell_end(&self, i: usize) -> usize {
        let mut end = i;
        while self.nest.ptn[end] > self.level {
            end += 1;
        }
        end
    }

    pub fn raw_cells(&self) -> PartitionNestChunkBy<'_> {
        self.nest.chunk_by(self.level)
    }

    pub fn cells(&self) -> PartitionCellsIter<'_> {
        PartitionCellsIter {
            partition: self,
            index: 0,
        }
    }

    pub fn cells_mut(&mut self) -> PartitionCellsIterMut<'_> {
        PartitionCellsIterMut {
            partition: self,
            index: 0,
        }
    }

    pub fn split(&mut self, i: usize) {
        self.nest.lab[i] = self.level;
        self.nest.numcells[self.level] += 1;
    }

    pub(crate) fn len(&self) -> usize {
        self.nest.lab.len()
    }
}

impl Index<usize> for Partition {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.nest.lab[index]
    }
}

pub struct Cell<'a> {
    pub(crate) partition: &'a Partition,
    pub(crate) first_lab_index: usize,
    pub(crate) cell_lab: &'a [usize],
}

impl<'a> Cell<'a> {
    pub fn is_at_end(&self) -> bool {
        self.first_lab_index == self.partition.nest.lab.len()
    }

    pub fn len(&self) -> usize {
        self.cell_lab.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.cell_lab.is_empty()
    }
}

/// https://doc.rust-lang.org/src/core/iter/traits/iterator.rs.html#2334-2337
///
/// Reorders the elements of this iterator *in-place* according to the given predicate,
/// such that all those that return `true` precede all those that return `false`.
/// Returns the number of `true` elements found.
///
/// The relative order of partitioned items is not maintained.
///
/// # Current implementation
///
/// The current algorithm tries to find the first element for which the predicate evaluates
/// to false and the last element for which it evaluates to true, and repeatedly swaps them.
///
/// Time complexity: *O*(*n*)
///
/// See also [`is_partitioned()`] and [`partition()`].
///
/// [`is_partitioned()`]: Iterator::is_partitioned
/// [`partition()`]: Iterator::partition
///
/// # Examples
///
/// ```
/// #![feature(iter_partition_in_place)]
///
/// let mut a = [1, 2, 3, 4, 5, 6, 7];
///
/// // Partition in-place between evens and odds
/// let i = a.iter_mut().partition_in_place(|n| n % 2 == 0);
///
/// assert_eq!(i, 3);
/// assert!(a[..i].iter().all(|n| n % 2 == 0)); // evens
/// assert!(a[i..].iter().all(|n| n % 2 == 1)); // odds
/// ```
pub(crate) fn partition_in_place<P>(slice: &mut [usize], ref mut predicate: P) -> usize
where
    P: FnMut(&usize) -> bool,
{
    // FIXME: should we worry about the count overflowing? The only way to have more than
    // `usize::MAX` mutable references is with ZSTs, which aren't useful to partition...

    // These closure "factory" functions exist to avoid genericity in `Self`.

    #[inline]
    fn is_false<'a, T>(
        predicate: &'a mut impl FnMut(&T) -> bool,
        true_count: &'a mut usize,
    ) -> impl FnMut(&&mut T) -> bool + 'a {
        move |x| {
            let p = predicate(&**x);
            *true_count += p as usize;
            !p
        }
    }

    #[inline]
    fn is_true<T>(predicate: &mut impl FnMut(&T) -> bool) -> impl FnMut(&&mut T) -> bool + '_ {
        move |x| predicate(&**x)
    }

    let mut iter = slice.into_iter();
    // Repeatedly find the first `false` and swap it with the last `true`.
    let mut true_count = 0;
    while let Some(head) = iter.find(is_false(predicate, &mut true_count)) {
        if let Some(tail) = iter.rfind(is_true(predicate)) {
            std::mem::swap(head, tail);
            true_count += 1;
        } else {
            break;
        }
    }
    true_count
}

pub struct CellMut<'a> {
    pub(crate) level: usize,
    pub(crate) first_lab_index: usize,
    pub(crate) cell_lab: &'a mut [usize],
    pub(crate) cell_ptn: &'a mut [usize],
    pub(crate) numcells: &'a mut usize,
}

impl<'a> CellMut<'a> {
    pub fn len(&self) -> usize {
        self.cell_lab.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.cell_lab.is_empty()
    }

    pub fn swap(&mut self, a: usize, b: usize) {
        self.cell_lab.swap(a, b)
    }

    pub fn split(&mut self, i: usize) {
        self.cell_ptn[i] = self.level;
        *self.numcells += 1;
    }

    pub fn partition_index(&self, i: usize) -> usize {
        self.first_lab_index + i
    }

    /// partition the cell, put splitter neighbours first
    /// returns the number of neighbours
    pub fn partition_in_place_std(&mut self, splitter_neighbours: &Set) -> usize {
        partition_in_place(self.cell_lab, |c| splitter_neighbours[*c])
    }

    /// partition the cell, put splitter neighbours first
    /// returns the number of neighbours
    pub fn partition_in_place_orig(&mut self, splitter_neighbours: &Set) -> usize {
        let mut c1: usize = 0;
        // c2 can be -1
        let mut c2: isize = self.len() as isize - 1;
        while c1 as isize <= c2 {
            if splitter_neighbours[self[c1]] {
                c1 += 1;
            } else {
                self.swap(c1, c2 as usize);
                c2 -= 1;
            }
        }
        c1
    }
}

impl<'a> Index<usize> for CellMut<'a> {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.cell_lab[index]
    }
}

impl<'a> IndexMut<usize> for CellMut<'a> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.cell_lab[index]
    }
}

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionCellsIter<'a> {
    pub(crate) partition: &'a Partition,
    pub(crate) index: usize,
}

impl<'a> Iterator for PartitionCellsIter<'a> {
    type Item = Cell<'a>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.index == self.partition.nest.lab.len() {
            None
        } else {
            let mut len = 1;
            let ptn_iter = self.partition.nest.ptn[self.index..].iter();
            for ptn in ptn_iter {
                if *ptn > self.partition.level {
                    len += 1
                } else {
                    break;
                }
            }
            let first_lab_index = self.index;
            self.index += len;
            Some(Cell {
                partition: self.partition,
                first_lab_index,
                cell_lab: &self.partition.nest.lab[first_lab_index..self.index],
            })
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.index == self.partition.nest.lab.len() {
            (0, Some(0))
        } else {
            (1, Some(self.partition.nest.lab.len() - self.index))
        }
    }

    #[inline]
    fn last(self) -> Option<Self::Item> {
        todo!()
    }
}

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionCellsIterMut<'a> {
    pub(crate) partition: &'a mut Partition,
    pub(crate) index: usize,
}

// do not implement trait Iterator because of lifetime issue
impl<'a> PartitionCellsIterMut<'a> {
    pub fn next<'b>(&'b mut self) -> Option<CellMut<'b>> {
        if self.index == self.partition.len() {
            None
        } else {
            let mut len = 1;
            for ptn in self.partition.nest.ptn[self.index..].iter() {
                if *ptn > self.partition.level {
                    len += 1
                } else {
                    break;
                }
            }
            let first_lab_index = self.index;
            self.index += len;
            let lab_slice: &'b mut [usize] =
                &mut self.partition.nest.lab[first_lab_index..self.index];
            let ptn_slice: &'b mut [usize] =
                &mut self.partition.nest.ptn[first_lab_index..self.index];
            Some(CellMut {
                level: self.partition.level,
                first_lab_index,
                cell_lab: lab_slice,
                cell_ptn: ptn_slice,
                numcells: &mut self.partition.nest.numcells[self.partition.level],
            })
        }
    }
}
