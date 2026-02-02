use crate::{graph6::BIT, gtools::g6string::G6String};

/// From nauty.h:
///
/// *   Conventions and Assumptions:                                             *
/// *                                                                            *
/// *    A 'setword' is the chunk of memory that is occupied by one part of      *
/// *    a set.  This is assumed to be >= WORDSIZE bits in size.                 *
/// *                                                                            *
/// *    The rightmost (loworder) WORDSIZE bits of setwords are numbered         *
/// *    0..WORDSIZE-1, left to right.  It is necessary that the 2^WORDSIZE      *
/// *    setwords with the other bits zero are totally ordered under <,=,>.      *
/// *    This needs care on a 1's-complement machine.                            *
/// *                                                                            *
/// *    The int variables m and n have consistent meanings throughout.          *
/// *    Graphs have n vertices always, and sets have m setwords always.         *
/// *                                                                            *
/// *    A 'set' consists of m contiguous setwords, whose bits are numbered      *
/// *    0,1,2,... from left (high-order) to right (low-order), using only       *
/// *    the rightmost WORDSIZE bits of each setword.  It is used to             *
/// *    represent a subset of {0,1,...,n-1} in the usual way - bit number x     *
/// *    is 1 iff x is in the subset.  Bits numbered n or greater, and           *
/// *    unnumbered bits, are assumed permanently zero.                          *
/// *                                                                            *
/// *    A 'graph' consists of n contiguous sets.  The i-th set represents       *
/// *    the vertices adjacent to vertex i, for i = 0,1,...,n-1.                 *
/// *                                                                            *
/// *    A 'permutation' is an array of n ints repesenting a permutation of      *
/// *    the set {0,1,...,n-1}.  The value of the i-th entry is the number to    *
/// *    which i is mapped.                                                      *
/// *                                                                            *
/// *    If g is a graph and p is a permutation, then g^p is the graph in        *
/// *    which vertex i is adjacent to vertex j iff vertex p[i] is adjacent      *
/// *    to vertex p[j] in g.                                                    *
/// *                                                                            *
/// *    A partition nest is represented by a pair (lab,ptn), where lab and ptn  *
/// *    are int arrays.  The "partition at level x" is the partition whose      *
/// *    cells are {lab[i],lab[i+1],...,lab[j]}, where [i,j] is a maximal        *
/// *    subinterval of [0,n-1] such that ptn[k] > x for i <= k < j and          *
/// *    ptn[j] <= x.  The partition at level 0 is given to nauty by the user.   *
/// *    This is  refined for the root of the tree, which has level 1.           *
///
/// here WORDSIZE == 64
const WORDSIZE: usize = 64;
const LOG_WORDSIZE: u8 = 6;
pub struct Graph<const M: usize>(pub Vec<Set<M>>);

/// the vertices adjacent to vertex i
pub struct Set<const M: usize> {
    pub words: [SetWord; M],
}

pub type SetWord = u64;

impl<const M: usize> Graph<M> {
    pub fn n(&self) -> usize {
        self.0.len()
    }

    // function ntog6 in gtools.c in nauty
    // https://users.cecs.anu.edu.au/~bdm/data/formats.txt
    pub fn to_graph6(&self) -> String {
        let mut g6 = G6String::from(self);
        g6.finish();
        g6.to_string()
    }
}

impl<const M: usize> Set<M> {
    pub fn is_element(&self, pos: usize) -> bool {
        (self.words[Self::set_wd(pos)] & BIT[Self::set_bt(pos)]) != 0
    }

    pub fn set_wd(pos: usize) -> usize {
        pos >> LOG_WORDSIZE
    }

    pub fn set_bt(pos: usize) -> usize {
        pos & (WORDSIZE - 1)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bit() {
        // for i in 0.. 10000 {
        //     let j =
        // }
        assert_eq!(64, BIT.len());
        for (i_bit, bit) in BIT.into_iter().enumerate() {
            assert_eq!(1 << (63 - i_bit), bit)
        }
    }
    // example from https://users.cecs.anu.edu.au/~bdm/data/formats.txt, line 73
    #[test]
    fn test_to_graph6() {
        let g = create_example();
        let g6 = g.to_graph6();
        assert_eq!(g6.bytes().collect::<Vec<_>>(), [68, 81, 99, 10]);
        assert_eq!(g6, "DQc\n");
    }

    #[test]
    fn test_is_element() {
        let g = create_example();
        assert!(g.0[3].is_element(1));
    }

    // example from https://users.cecs.anu.edu.au/~bdm/data/formats.txt, line 73
    fn create_example() -> Graph<1> {
        Graph(vec![
            Set {
                words: [BIT[2] | BIT[4]], // 0 0-2, 0-4
            },
            Set {
                words: [BIT[3]], // 1 1-3
            },
            Set {
                words: [BIT[0]], // 2 0-2
            },
            Set {
                words: [BIT[1] | BIT[4]], // 3  1-3 3-4
            },
            Set {
                words: [BIT[0] | BIT[3]], // 4 0-4 3-4
            },
        ])
    }
}
