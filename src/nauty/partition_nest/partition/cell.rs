use super::Partition;
use crate::nauty::{Graph, Set, SetTrait};
use bitvec::{bitvec, order::Msb0};
use std::ops::Index;

#[derive(Clone)]
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

    pub fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.cell_lab.iter()
    }

    pub fn to_vec(&self) -> Vec<usize> {
        self.cell_lab.to_vec()
    }

    pub fn get_splitters(&self) -> Set {
        let mut set = bitvec![usize, Msb0; 0; self.partition.nest.lab.len()];
        for i in self.cell_lab {
            set.set(*i, true);
        }
        set
    }

    pub fn set(&self, n: usize) -> Set {
        let mut cell_set = Set::zeros(n);
        for i in self.iter() {
            cell_set.add_one(*i);
        }
        cell_set
    }
}

impl<'a> Index<usize> for Cell<'a> {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.cell_lab[index]
    }
}
