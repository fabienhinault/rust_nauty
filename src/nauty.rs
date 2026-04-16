use std::{
    collections::HashSet,
    ops::{Index, IndexMut},
};

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
pub type SetWord = BitVec<usize, Msb0>;
#[derive(Default)]
pub struct Graph(pub Vec<BitVec<usize, Msb0>>);
pub type NautyCounter = u128;

trait SetWordTrait {
    fn difference(&self, other: &Self) -> Self;
    fn first_bit_nz_index(&self) -> usize;
    fn one(index: usize) -> Self;
    fn add_one(&mut self, index: usize);
    fn remove_one(&mut self, index: usize);
    fn except_one(&self, index: usize) -> Self;
    fn filter(&self, other: &Self) -> Self;
}

impl SetWordTrait for SetWord {
    fn difference(&self, other: &Self) -> Self {
        self.clone() & !other.clone()
    }

    fn add_one(&mut self, index: usize) {
        *self |= Self::one(index);
    }

    fn remove_one(&mut self, index: usize) {
        *self &= !Self::one(index);
    }

    fn except_one(&self, index: usize) -> Self {
        self.difference(&Self::one(index))
    }

    fn first_bit_nz_index(&self) -> usize {
        self.leading_zeros()
    }

    fn one(index: usize) -> Self {
        Self::from_element(BIT[index])
    }

    fn filter(&self, other: &Self) -> Self {
        self.clone() & other
    }
}

impl IndexMut<usize> for Graph {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl Index<usize> for Graph {
    type Output = SetWord;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Graph {
    // graph with one vertex and no edge
    pub fn one() -> Self {
        Self(vec![bitvec![usize, Msb0; 0]])
    }

    pub fn n(&self) -> usize {
        self.0.len()
    }

    pub fn get(&self, i: usize) -> &SetWord {
        &self.0[i]
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
    // The algorithm a subpart of the finding of biconnected components in https://dl.acm.org/doi/epdf/10.1145/362248.362272
    //
    // The algorithm terminates in 2 * V max.
    // Proof: It passes no more than V times in the 'if' because each time, a vertex is added to visited, and when all vertices
    // are in visited, it will not pass in the 'if' any more.
    // It passes no more than V times in the 'else', because each time a vertex is popped from it if it does not return,
    // and the stack size grows to V at maximum.
    //
    // The algorithm is correct, i.e. it returns true iff the graph is biconnected.
    /* test if g is biconnected */
    pub fn isbiconnected(&self) -> bool {
        let n = self.n();
        if n <= 2 {
            return false;
        }
        let mut visited = BitVec::from_element(BIT[0]);
        let mut stack = vec![0];
        // Numbering of vertices by order of discovery
        let mut discovery = VecMap(vec![Some(0)]);
        // For each vertex, the lowest point on the stack to which it is connected by another path of visited points (not the one of the DFS)
        // This value is set progressively.
        // Only for neighbours of the lowest point during the pass in 'if', then during the pass in 'else' for the other connected points.
        // The value is the good one just at time for the test 'return false'.
        let mut low_point = VecMap(vec![Some(0)]);
        let mut numvis = 1_usize;
        let mut v = 0_usize;
        let mut w;

        loop {
            let not_visited = self[v].difference(&visited);
            if not_visited.any() {
                w = v;
                v = not_visited.first_bit_nz_index(); /* visit next child */
                stack.push(v);
                visited.add_one(v);
                numvis += 1;
                low_point.set(v, numvis);
                discovery.set(v, numvis);
                let mut visited_adjascents_not_parent = self[v].filter(&visited).except_one(w);
                while visited_adjascents_not_parent.any() {
                    w = visited_adjascents_not_parent.first_bit_nz_index();
                    visited_adjascents_not_parent.remove_one(w);
                    if discovery.get(w) < low_point.get(v) {
                        low_point.set(v, discovery.get(w));
                    }
                }
            } else {
                w = v; /* back up to parent */
                if stack.len() <= 1 {
                    // Visited the whole connected component containing 0, found no articulation point.
                    // biconnected iff visited whole graph.
                    return numvis == n;
                }
                v = stack.pop().unwrap();
                if low_point.get(w) >= discovery.get(v) {
                    return false;
                }
                if low_point.get(w) < low_point.get(v) {
                    low_point.set(v, low_point.get(w));
                }
            }
        }
    }
}

struct VecMap(Vec<Option<usize>>);

impl VecMap {
    pub fn new() -> Self {
        Self(vec![])
    }

    pub fn get(&self, index: usize) -> usize {
        self.0[index].unwrap()
    }

    pub fn set(&mut self, index: usize, value: usize) {
        while self.0.len() < index + 1 {
            self.0.push(None);
        }
        self.0[index] = Some(value);
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use test_case::test_case;

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

    #[test_case(create_n_circle(6))]
    #[test_case(create_diamond())]
    #[test_case(create_complete(4))]
    #[test_case(create_g4g_bc())]
    fn test_is_biconnected(biconnected_graph: Graph) {
        assert!(biconnected_graph.isbiconnected());
    }

    #[test_case(create_n_path(4))]
    #[test_case(create_g4g_not_bc())]
    #[test_case(create_a())]
    fn test_is_not_biconnected(not_biconnected_graph: Graph) {
        assert!(!not_biconnected_graph.isbiconnected());
    }

    #[test]
    fn test_diamond_isbiconnected() {
        assert!(create_diamond().isbiconnected());
    }

    #[test]
    fn test_tetraedron_isbiconnected() {
        assert!(create_complete(4).isbiconnected());
    }

    #[test]
    fn test_g4g_not_isbiconnected() {
        assert!(!create_g4g_not_bc().isbiconnected());
    }

    #[test]
    fn test_without_loop_not_isbiconnected() {
        assert!(!create_without_loop().isbiconnected());
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

    //  0---1---2---3---4--- ... ---n
    fn create_n_path(n: usize) -> Graph {
        let mut first = bitvec![usize, Msb0; 0; n];
        first.set(1, true);
        let mut last = bitvec![usize, Msb0; 0; n];
        last.set(n - 2, true);
        let mut words = vec![first];
        words.extend((1..n - 1).map(|i| create_n_circle_bitvec(n, i)));
        words.push(last);
        Graph(words)
    }

    fn create_diamond() -> Graph {
        Graph(vec![
            //                   0  1  2  3
            bitvec![usize, Msb0; 0, 1, 1, 1],
            bitvec![usize, Msb0; 1, 0, 1, 0],
            bitvec![usize, Msb0; 1, 1, 0, 1],
            bitvec![usize, Msb0; 1, 0, 1, 0],
        ])
    }

    fn create_complete(n: usize) -> Graph {
        Graph((0..n).map(|i| create_complete_bitvec(n, i)).collect())
    }

    fn create_complete_bitvec(n: usize, i: usize) -> SetWord {
        let mut result: BitVec<usize, Msb0> = bitvec![usize, Msb0; 1; n];
        result.set(i, false);
        result
    }

    // https://www.geeksforgeeks.org/dsa/biconnectivity-in-a-graph/
    // 1-0--3
    // |/   |
    // 2    4
    fn create_g4g_not_bc() -> Graph {
        Graph(vec![
            //                   0  1  2  3  4
            bitvec![usize, Msb0; 0, 1, 1, 1, 0],
            bitvec![usize, Msb0; 1, 0, 1, 0, 0],
            bitvec![usize, Msb0; 1, 1, 0, 0, 0],
            bitvec![usize, Msb0; 1, 0, 0, 0, 1],
            bitvec![usize, Msb0; 0, 0, 0, 1, 0],
        ])
    }

    //  1-0--3
    //  |/   |
    //  2----4
    fn create_g4g_bc() -> Graph {
        Graph(vec![
            //                   0  1  2  3  4
            bitvec![usize, Msb0; 0, 1, 1, 1, 0],
            bitvec![usize, Msb0; 1, 0, 1, 0, 0],
            bitvec![usize, Msb0; 1, 1, 0, 0, 1],
            bitvec![usize, Msb0; 1, 0, 0, 0, 1],
            bitvec![usize, Msb0; 0, 0, 1, 1, 0],
        ])
    }

    //  ,----.
    //  1-0--3
    //  |/   |
    //  2    4
    fn create_a() -> Graph {
        Graph(vec![
            //                   0  1  2  3  4
            bitvec![usize, Msb0; 0, 1, 1, 1, 0],
            bitvec![usize, Msb0; 1, 0, 1, 1, 0],
            bitvec![usize, Msb0; 1, 1, 0, 0, 0],
            bitvec![usize, Msb0; 1, 1, 0, 0, 1],
            bitvec![usize, Msb0; 0, 0, 0, 1, 0],
        ])
    }

    fn create_without_loop() -> Graph {
        Graph(vec![
            //                   0  1  2  3
            bitvec![usize, Msb0; 0, 1, 1, 1],
            bitvec![usize, Msb0; 1, 0, 0, 0],
            bitvec![usize, Msb0; 1, 0, 0, 0],
            bitvec![usize, Msb0; 1, 0, 0, 0],
        ])
    }
}
