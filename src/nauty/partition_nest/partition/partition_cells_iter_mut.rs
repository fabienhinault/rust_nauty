use std::num::NonZero;

use super::{Partition, cell_mut::CellMut};

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionCellsIterMut<'a> {
    pub(crate) partition: &'a mut Partition,
    pub(crate) index: usize,
}

// do not implement trait Iterator because of lifetime issue
// between self.partition's lifetime and self in next().
// ChunkByMut https://doc.rust-lang.org/stable/src/core/slice/iter.rs.html#3120
// gets out of it by splitting self.slice
// we could do the same for lab and ptn, like in PartitionCellsIter.
// But numcells remains. (put it in a RefCell? Or remove it from PartitionNest and compute it when needed?)
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

    pub fn nth<'b>(&'b mut self, n: usize) -> Option<CellMut<'b>> {
        for _i in 0..n {
            self.next()?;
        }
        self.next()
    }
}
