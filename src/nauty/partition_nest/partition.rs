use super::PartitionNest;
use super::partition_nest_chunk_by::PartitionNestChunkBy;
use crate::nauty::Set;
use crate::nauty::partition_nest::partition::cell::Cell;
use crate::nauty::partition_nest::partition::cell_mut::CellMut;
use bitvec::bitvec;
use bitvec::order::Msb0;
use std::ops::Index;

mod cell;
pub mod cell_mut;
mod partition_cells_iter;
mod partition_cells_iter_mut;
mod partition_subset;

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

    pub fn string(&self) -> String {
        self.nest.partition_string(self.level)
    }

    pub fn numcells(&self) -> usize {
        self.nest.numcells(self.level)
    }

    pub fn numcells_mut(&mut self) -> &mut usize {
        self.nest.numcells_mut(self.level)
    }

    pub fn is_discrete(&self) -> bool {
        self.nest.ptn.iter().all(|i| *i <= self.level)
    }

    /// last index of cell containing i for partition of given level
    pub fn cell_end(&self, i: usize) -> usize {
        let mut end = i;
        while self.nest.ptn[end] > self.level {
            end += 1;
        }
        end
    }

    pub fn get_splitters(&self, start: usize) -> Set {
        let mut set = bitvec![usize, Msb0; 0; self.nest.lab.len()];
        for i in start..=self.cell_end(start) {
            set.set(self.nest.lab[i], true);
        }
        set
    }

    pub fn raw_cells(&self) -> PartitionNestChunkBy<'_> {
        self.nest.chunk_by(self.level)
    }

    pub fn cells(&self) -> partition_cells_iter::PartitionCellsIter<'_> {
        partition_cells_iter::PartitionCellsIter {
            partition: self,
            index: 0,
        }
    }

    pub fn cells_mut(&mut self) -> partition_cells_iter_mut::PartitionCellsIterMut<'_> {
        partition_cells_iter_mut::PartitionCellsIterMut {
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

    pub fn get_cell<'a>(&'a self, i: usize) -> Cell<'a> {
        assert!(i == 0 || self.nest.ptn[i - 1] >= self.level);
        Cell {
            partition: self,
            first_lab_index: i,
            cell_lab: &self.nest.lab[i..=self.cell_end(i)],
        }
    }

    pub fn get_cell_mut<'a>(&'a mut self, i: usize) -> CellMut<'a> {
        assert!(i == 0 || self.nest.ptn[i - 1] >= self.level);
        let cell_end = self.cell_end(i);
        CellMut {
            level: self.level,
            first_lab_index: i,
            cell_lab: &mut self.nest.lab[i..=cell_end],
            cell_ptn: &mut self.nest.ptn[i..=cell_end],
            numcells: &mut self.nest.numcells[self.level],
        }
    }
}

impl Index<usize> for Partition {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.nest.lab[index]
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::nauty::partition_nest::partition::Partition;

    const LAB: &[usize] = &[4, 6, 2, 0, 8, 7, 5, 9, 3, 1];
    const PTN: &[usize] = &[4, 3, 4, 1, 2, 4, 4, 0, 2, 0];

    #[test]
    fn test_is_discrete_4() {
        let partition4 = Partition::new(PartitionNest::new(LAB.to_vec(), PTN.to_vec()), 4);
        assert!(partition4.is_discrete());
        assert_eq!(partition4.numcells(), partition4.len());
        assert_eq!(partition4.string(), "4 | 6 | 2 | 0 | 8 | 7 | 5 | 9 | 3 | 1");
    }

    #[test]
    fn test_is_not_discrete_3() {
        let partition3 = Partition::new(PartitionNest::new(LAB.to_vec(), PTN.to_vec()), 3);
        assert!(!partition3.is_discrete());
        assert_eq!(partition3.numcells(), 6);
        assert_eq!(partition3.string(), "4, 6 | 2, 0 | 8 | 7, 5, 9 | 3 | 1");
    }
}
