use super::{Partition, cell::Cell};

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
            for ptn in self.partition.nest.ptn[self.index..].iter() {
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
