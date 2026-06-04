use bitvec::{bitvec, order::Msb0};

use crate::nauty::{Set, SetTrait};

pub trait SetWordNautilTrait {
    fn next_element(&self, pos: isize) -> isize;
}

impl SetWordNautilTrait for Set {
    /*****************************************************************************
     *                                                                            *
     *  nextelement(set1,m,pos) = the position of the first element in set set1   *
     *  which occupies a position greater than pos.  If no such element exists,   *
     *  the value is -1.  pos can have any value less than n, including negative  *
     *  values.                                                                   *
     *                                                                            *
     *  GLOBALS ACCESSED: none                                                    *
     *                                                                            *
     *****************************************************************************/
    fn next_element(&self, pos: isize) -> isize {
        let setwd: Set = self.filter(&self.bit_mask(pos.max(0) as usize));
        let lz = setwd.leading_zeros();
        if lz == setwd.len() { -1 } else { lz as isize }
    }
}
