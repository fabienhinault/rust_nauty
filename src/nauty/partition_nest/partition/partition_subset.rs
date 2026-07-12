use crate::{
    nautil::SetWordNautilTrait,
    nauty::{
        Set,
        partition_nest::partition::{Partition, cell::Cell},
    },
};

pub struct PartitionSubset<'a> {
    partition: &'a Partition,
    cells_indices: Set,
    hint: usize,
}

pub struct PartitionSubsetIterator<'a> {
    subset: PartitionSubset<'a>,
    current_cell_index: usize,
}

impl<'a> Iterator for PartitionSubsetIterator<'a> {
    type Item = Cell<'a>;

    /// rotate over cells in subset until subset is empty
    fn next(&mut self) -> Option<Self::Item> {
        if self.subset.cells_indices[self.current_cell_index] {
            Some(self.subset.partition.get_cell(self.current_cell_index))
        } else {
            let next_element = self
                .subset
                .cells_indices
                .next_element(Some(self.current_cell_index));
            self.current_cell_index = match next_element {
                Some(next_element) => next_element,
                None => match self.subset.cells_indices.first_one() {
                    Some(first_one) => first_one,
                    None => {
                        return None;
                    }
                },
            };
            Some(self.subset.partition.get_cell(self.current_cell_index))
        }
    }
}


