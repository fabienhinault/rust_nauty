use crate::{
    graph6::BIT,
    gtools::{
        g6char::G6Char,
        g6error::G6Error,
        g6string::{G6String, graph_size},
    },
    nauty::{
        Set, SetTrait, VecMap, bitvec_from_closure, bitvec_from_f,
        partition_nest::partition::Partition, u32_to_bitvec,
    },
};
use bitvec::{bitvec, order::Msb0, vec::BitVec};

#[derive(Default, PartialEq, Debug, Clone)]
pub struct Graph(pub Vec<BitVec<usize, Msb0>>);

impl Graph {
    // graph with one vertex and no edge
    pub fn one() -> Self {
        Self(vec![bitvec![usize, Msb0; 0; 1]])
    }

    pub fn new(vertex_count: usize) -> Self {
        Self::no_edge(vertex_count)
    }

    pub fn no_edge(vertex_count: usize) -> Self {
        Self(vec![bitvec![usize, Msb0; 0; vertex_count]; vertex_count])
    }

    pub fn from_f(vertex_count: usize, f: fn(usize, usize, usize) -> bool) -> Self {
        Graph(
            (0..vertex_count)
                .map(|i| bitvec_from_f(vertex_count, i, f))
                .collect(),
        )
    }

    pub fn from_edges(vertex_count: usize, edges: &[(usize, usize)]) -> Self {
        let mut g = Self::no_edge(vertex_count);
        for (v0, v1) in edges {
            g.0[*v0].add_one(*v1);
            g.0[*v1].add_one(*v0);
        }
        g
    }

    pub fn from_closure<F: Fn(usize, usize, usize) -> bool + Copy>(
        vertex_count: usize,
        f: F,
    ) -> Self {
        Graph(
            (0..vertex_count)
                .map(|i| bitvec_from_closure(vertex_count, i, f))
                .collect(),
        )
    }

    pub fn join(first: &Self, second: &Self) -> Self {
        let mut upper = first
            .0
            .clone()
            .into_iter()
            .map(|mut bv| {
                bv.extend(bitvec![usize, Msb0; 1; second.n()]);
                bv
            })
            .collect::<Vec<_>>();
        upper.extend(second.0.clone().into_iter().map(|bv| {
            let mut bv2 = bitvec![usize, Msb0; 1; first.n()];
            bv2.extend(bv.iter());
            bv2
        }));
        Graph(upper)
    }

    pub fn n(&self) -> usize {
        self.0.len()
    }

    pub fn get(&self, i: usize) -> &Set {
        &self.0[i]
    }

    pub fn upper(&self, i: usize) -> Set {
        self[i].masked(i)
    }

    // function ntog6 in gtools.c in nauty
    // https://users.cecs.anu.edu.au/~bdm/data/formats.txt
    pub fn to_graph6(&self) -> String {
        G6String::from(self).to_string()
    }

    pub fn to_matrix(&self) -> String {
        itertools::Itertools::intersperse(
            self.0.iter().map(|row| {
                let s: String = itertools::Itertools::intersperse(
                    row.iter()
                        .map(|other_vertex| if *other_vertex { "1" } else { "0" }),
                    " ",
                )
                .collect();
                s
            }),
            "\n".to_owned(),
        )
        .collect()
    }

    pub(crate) fn from_graph6(g6: String) -> Result<Self, G6Error> {
        let mut iter: std::slice::Iter<'_, u8> = g6.as_bytes().iter();
        let n = graph_size(&mut iter)?;
        let mut g = Self::no_edge(n);
        let mut g6_char = G6Char::try_from(*iter.next().ok_or(G6Error())?)?;
        for i_row in 1..n {
            for i_other_vertex in 0..i_row {
                if g6_char.is_empty() {
                    g6_char = G6Char::try_from(*iter.next().ok_or(G6Error())?)?;
                }
                if g6_char.pop_front() {
                    g.0[i_row].set(i_other_vertex, true);
                    g.0[i_other_vertex].set(i_row, true);
                }
            }
        }
        Ok(g)
    }

    pub(crate) fn from_u32(input: &[u32]) -> Self {
        let n = input.len();
        Graph(input.iter().map(|&u| u32_to_bitvec(u, n)).collect())
    }

    pub fn canonise(&self) -> Self {
        todo!();
        self.clone()
    }

    /*****************************************************************************
     *                                                                            *
     *  testcanlab(g,canong,lab,samerows,m,n) compares g^lab to canong,           *
     *  using an ordering which is immaterial since it's only used here.  The     *
     *  value returned is -1,0,1 if g^lab <,=,> canong.  *samerows is set to      *
     *  the number of rows (0..n) of canong which are the same as those of g^lab. *
     *                                                                            *
     *  GLOBALS ACCESSED: workset<rw>,permset(),workperm<rw>                      *
     *                                                                            *
     *****************************************************************************/

    pub fn testcanlab(&self, canong: &Graph, partition: &Partition, samerows: &mut usize) -> isize {
        let workperm = partition.permutation1();
        for i in 0..self.n() {
            let workset = self.0[partition[i]].permset(&workperm);
            if workset < canong.0[i] {
                *samerows = i;
                return -1;
            }
            if workset > canong.0[i] {
                *samerows = i;
                return 1;
            }
        }
        *samerows = self.n();
        return 0;
    }

    /*****************************************************************************
     *                                                                            *
     *  updatecan(g,canong,lab,samerows,m,n) sets canong = g^lab, assuming        *
     *  the first samerows of canong are ok already.                              *
     *                                                                            *
     *  GLOBALS ACCESSED: permset(),workperm<rw>                                  *
     *                                                                            *
     *****************************************************************************/
    // naugraph.c 143
    pub fn updatecan(&self, partition: &Partition, samerows: usize) -> Self {
        let workperm = partition.permutation1();
        let mut canon = self.clone();
        for i in samerows..self.n() {
            canon.0[i] = self.0[partition[i]].permset(&workperm);
        }
        canon
    }

    /* test if g is connected */
    // geng.c 636
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
