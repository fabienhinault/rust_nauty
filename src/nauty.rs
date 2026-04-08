use crate::{MAXN, graph6::BIT, gtools::g6string::G6String};
use bitvec::{bitvec, order::Msb0, vec::BitVec};

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
/// here WORDSIZE == size_of(usize) (64)
pub const WORDSIZE: usize = usize::BITS as usize;
const LOG_WORDSIZE: u8 = (WORDSIZE - 1).count_ones() as u8;

// the BitVec of index i has the vertices adjascent to vertex of index i.
// g.0[i][j] == 1 iff (i, j) is an edge of g.
#[derive(Default)]
pub struct Graph(pub Vec<BitVec<usize, Msb0>>);
pub type NautyCounter = u128;

impl Graph {
    // graph with one vertex and no edge
    pub fn one() -> Self {
        Self(vec![bitvec![usize, Msb0; 0]])
    }

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

    pub fn isconnected(&self) -> bool {
        let n = self.n();
        let allbits = bitvec![1; n];
        let mut expanded = BitVec::from_element(BIT[n - 1]);
        let mut seen = expanded.clone() | self.0[n - 1].clone();
        let mut toexpand = seen.clone() & !expanded.clone();
        while seen != allbits && toexpand.any() {
            let i = toexpand.leading_zeros();
            expanded |= BitVec::from_element(BIT[i]);
            seen |= self.0[i].clone();
            toexpand = seen.clone() & !expanded.clone();
        }
        seen[..n] == allbits
    }

    // static boolean isbiconnected(graph *g, int n)
    // https://en.wikipedia.org/wiki/Biconnected_graph
    // A connected graph that is not broken into disconnected pieces by deleting any single vertex (and incident edges).
    /* test if g is biconnected */
    pub fn isbiconnected(&self) -> bool {
        let n = self.n();
        if n <= 2 {
            return false;
        }
        let mut visited = BitVec::from_element(BIT[0]);
        let mut stack = vec![];
        let mut num = vec![0_usize; MAXN];
        let mut lp = vec![0_usize; MAXN];
        let mut numvis = 1_usize;
        let mut v = 0_usize;
        let mut w;

        loop {
            let mut sw = self.0[v].clone() & !visited.clone();
            if sw.any() {
                w = v;
                v = sw.leading_zeros(); /* visit next child */
                stack.push(v);
                visited |= BitVec::from_element(BIT[v]);
                numvis += 1;
                lp[v] = numvis;
                num[v] = numvis;
                sw = self.0[v].clone() & visited.clone() & !BitVec::from_element(BIT[w]);
                while sw.any() {
                    w = sw.leading_zeros();
                    sw &= !BitVec::from_element(BIT[w]);
                    if num[w] < lp[v] {
                        lp[v] = num[w];
                    }
                }
            } else {
                w = v; /* back up to parent */
                if stack.len() <= 1 {
                    return numvis == n;
                }
                v = stack.pop().unwrap();
                if lp[w] >= num[v] {
                    return false;
                }
                if lp[w] < lp[v] {
                    lp[v] = lp[w];
                }
            }
        }
    }
}

#[cfg(test)]
mod test {

    use super::*;

    // #[test]
    // fn test_bit() {
    //     // for i in 0.. 10000 {
    //     //     let j =
    //     // }
    //     assert_eq!(64, BIT.len());
    //     for (i_bit, bit) in BIT.into_iter().enumerate() {
    //         assert_eq!(1 << (63 - i_bit), bit)
    //     }
    // }
    // example from https://users.cecs.anu.edu.au/~bdm/data/formats.txt, line 73
    #[test]
    fn test_to_graph6() {
        let g = create_example();
        let g6 = g.to_graph6();
        assert_eq!(g6.bytes().collect::<Vec<_>>(), [68, 81, 99, 10]);
        assert_eq!(g6, "DQc\n");
    }

    #[test]
    fn test_example_is_connected() {
        assert!(create_example().isconnected());
    }

    #[test]
    fn test_disconnected() {
        assert!(!create_disconnected().isconnected());
    }

    #[test]
    fn test_disconnected_isbiconnected() {
        assert!(!create_disconnected().isbiconnected());
    }

    #[test]
    fn test_example_isbiconnected() {
        assert!(!create_example().isbiconnected());
    }

    #[test]
    fn test_circle_isbiconnected() {
        assert!(create_circle().isbiconnected());
    }

    #[test]
    fn test_n_circle_isbiconnected() {
        assert!(create_n_circle(6).isbiconnected());
    }

    // example from https://users.cecs.anu.edu.au/~bdm/data/formats.txt, line 73
    //
    //  2---0---4---3---1
    //
    fn create_example() -> Graph {
        Graph(vec![
            //                   0  1  2  3  4
            bitvec![usize, Msb0; 0, 0, 1, 0, 1],
            bitvec![usize, Msb0; 0, 0, 0, 1, 0],
            bitvec![usize, Msb0; 1, 0, 0, 0, 0],
            bitvec![usize, Msb0; 0, 1, 0, 0, 1],
            bitvec![usize, Msb0; 1, 0, 0, 1, 0],
        ])
    }

    //  2---0   4---3---1
    fn create_disconnected() -> Graph {
        Graph(vec![
            //                   0  1  2  3  4
            bitvec![usize, Msb0; 0, 0, 1, 0, 0],
            bitvec![usize, Msb0; 0, 0, 0, 1, 0],
            bitvec![usize, Msb0; 1, 0, 0, 0, 0],
            bitvec![usize, Msb0; 0, 1, 0, 0, 1],
            bitvec![usize, Msb0; 0, 0, 0, 1, 0],
        ])
    }

    //  ,---------------.
    //  2---0---4---3---1
    fn create_circle() -> Graph {
        Graph(vec![
            //                   0  1  2  3  4
            bitvec![usize, Msb0; 0, 0, 1, 0, 1],
            bitvec![usize, Msb0; 0, 0, 1, 1, 0],
            bitvec![usize, Msb0; 1, 1, 0, 0, 0],
            bitvec![usize, Msb0; 0, 1, 0, 0, 1],
            bitvec![usize, Msb0; 1, 0, 0, 1, 0],
        ])
    }

    //  ,---------------- ... --.
    //  0---1---2---3---4 ... --n
    fn create_n_circle(n: usize) -> Graph {
        Graph((0..n).map(|i| create_n_circle_bitvec(n, i)).collect())
    }

    fn create_n_circle_bitvec(n: usize, i: usize) -> BitVec<usize, Msb0> {
        let mut result: BitVec<usize, Msb0> = bitvec![usize, Msb0; 0; n];
        result.set((i + 1) % n, true);
        let i = i as isize;
        let n = n as isize;
        let im1 = (i - 1).rem_euclid(n);
        result.set(im1 as usize, true);
        result
    }
}
