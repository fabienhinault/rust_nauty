use super::Partition;

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
