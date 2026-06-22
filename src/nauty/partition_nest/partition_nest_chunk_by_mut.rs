use std::mem;

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionNestChunkByMut<'a> {
    lab_slice: &'a mut [usize],
    ptn_slice: &'a mut [usize],
    level: usize,
}

impl<'a> PartitionNestChunkByMut<'a> {
    pub fn new(lab_slice: &'a mut [usize], ptn_slice: &'a mut [usize], level: usize) -> Self {
        Self {
            lab_slice,
            ptn_slice,
            level,
        }
    }
}

impl<'a> Iterator for PartitionNestChunkByMut<'a> {
    type Item = &'a mut [usize];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.lab_slice.is_empty() {
            None
        } else {
            let mut len = 1;
            let ptn_iter = self.ptn_slice.iter();
            for ptn in ptn_iter {
                if *ptn > self.level {
                    len += 1
                } else {
                    break;
                }
            }
            let lab_slice: &'a mut [usize] = mem::take(&mut self.lab_slice);
            let ptn_slice: &'a mut [usize] = mem::take(&mut self.ptn_slice);
            let (lab_head, lab_tail) = lab_slice.split_at_mut(len);
            self.lab_slice = lab_tail;
            let (_ptn_head, ptn_tail) = ptn_slice.split_at_mut(len);
            self.ptn_slice = ptn_tail;
            Some(lab_head)
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.lab_slice.is_empty() {
            (0, Some(0))
        } else {
            (1, Some(self.lab_slice.len()))
        }
    }

    #[inline]
    fn last(mut self) -> Option<Self::Item> {
        todo!()
    }
}
