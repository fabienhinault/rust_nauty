use crate::nauty::VecMap;

use super::Partition;
use super::cell_mut::CellMut;
use std::ops::Index;

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

impl<'a> Index<usize> for Cell<'a> {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.cell_lab[index]
    }
}
