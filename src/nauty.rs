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
/// here WORDSIZE == size_of(usize)
use crate::{
    graph6::BIT,
    gtools::{
        g6char::G6Char,
        g6error::G6Error,
        g6string::{G6String, graph_size},
    },
};
use bitvec::{bitvec, order::Msb0, vec::BitVec, view::BitView};
use std::fmt::Debug;
use std::{
    collections::LinkedList,
    fs::File,
    ops::{Index, IndexMut},
};

pub mod partition_nest;
#[cfg(test)]
pub mod test;

struct OptionBlk {
    getcanon: u8,       /* make canong and canonlab? */
    digraph: bool,      /* multiple edges or loops? */
    writeautoms: bool,  /* write automorphisms? */
    writemarkers: bool, /* write stats on pts fixed, etc.? */
    defaultptn: bool,   /* set lab,ptn,active for single cell? */
    cartesian: bool,    /* use cartesian rep for writing automs? */
    linelength: u8,     /* max chars/line (excl. '\n') for output */
    outfile: File,      /* file for output, if any */
    tc_level: u8,       /* max level for smart target cell choosing */
    mininvarlevel: u8,  /* min level for invariant computation */
    maxinvarlevel: u8,  /* max level for invariant computation */
    invararg: u8,       /* value passed to (*invarproc)() */
    schreier: bool,     /* use random schreier method */  // skip for now
}

struct StatBlk {
    grpsize1: f64,        /* size of group is */
    grpsize2: i64,        /*    grpsize1 * 10^grpsize2 */
    numorbits: usize,     /* number of orbits in group */
    numgenerators: usize, /* number of generators found */
    errstatus: u8,        /* if non-zero : an error code */
    numnodes: usize,      /* total number of nodes */
    numbadleaves: usize,  /* number of leaves of no use */
    maxlevel: usize,      /* maximum depth of search */
    tctotal: usize,       /* total size of all target cells */
    canupdates: usize,    /* number of updates of best label */
    invapplics: usize,    /* number of applications of invarproc */
    invsuccesses: usize,  /* number of successful uses of invarproc() */
    invarsuclevel: usize, /* least level where invarproc worked */
}

impl StatBlk {
    fn new(numorbits: usize) -> Self {
        Self {
            numorbits,
            ..Default::default()
        }
    }
}

impl Default for StatBlk {
    fn default() -> Self {
        Self {
            grpsize1: 1.0,
            grpsize2: 0,
            numorbits: 0,
            numgenerators: 0,
            errstatus: 0,
            numnodes: 0,
            numbadleaves: 0,
            maxlevel: 0,
            tctotal: 0,
            canupdates: 0,
            invapplics: 0,
            invsuccesses: 0,
            invarsuclevel: 0,
        }
    }
}

pub const WORDSIZE: usize = usize::BITS as usize;
const LOG_WORDSIZE: u8 = (WORDSIZE - 1).count_ones() as u8;
pub const NAUTY_INFINITY: usize = 2_000_000_002; /* Max graph size is 2 billion */
pub const NAUTY_INFINITY_I: isize = 2_000_000_002; /* Max graph size is 2 billion */

// the BitVec of index i has the vertices adjascent to vertex of index i.
// g.0[i][j] == 1 iff (i, j) is an edge of g.
pub type Set = BitVec<usize, Msb0>;
#[derive(Default, PartialEq, Debug)]
pub struct Graph(pub Vec<BitVec<usize, Msb0>>);
pub type NautyCounter = u128;

pub trait SetTrait {
    fn difference(&self, other: &Self) -> Self;
    fn first_bit_nz_index(&self) -> usize;
    fn one(index: usize) -> Self;
    fn add_one(&mut self, index: usize);
    fn remove_one(&mut self, index: usize);
    fn except_one(&self, index: usize) -> Self;
    fn filter(&self, other: &Self) -> Self;
    fn ones_iter(&self) -> SetIterator;
    fn rounding_ones_iter(&self, hint: usize) -> RoundingSetIterator;
    fn bit_mask(&self, pos: usize) -> Self;
    fn masked(&self, pos: usize) -> Self;
}

impl SetTrait for Set {
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

    fn ones_iter(&self) -> SetIterator {
        SetIterator { set: self.clone() }
    }

    fn rounding_ones_iter(&self, hint: usize) -> RoundingSetIterator {
        RoundingSetIterator {
            set: self.clone(),
            current_index: hint,
        }
    }

    fn bit_mask(&self, pos: usize) -> Self {
        let mut result = bitvec![usize, Msb0; 0; pos];
        result.extend_from_bitslice(self[0..self.len() - pos].iter().as_bitslice());
        result
    }

    fn masked(&self, pos: usize) -> Self {
        self.filter(&self.bit_mask(pos))
    }
}

pub struct RoundingSetIterator {
    set: Set,
    current_index: usize,
}

impl Iterator for RoundingSetIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let setwd: Set = self.set.bit_mask(self.current_index);
        let mut lz = setwd.leading_zeros();
        if lz == setwd.len() {
            lz = self.set.leading_zeros();
        }
        if lz == self.set.len() {
            return None;
        }
        self.current_index = (lz + 1) % self.set.len();
        Some(lz)
    }
}

pub struct SetIterator {
    set: Set,
}

impl Iterator for SetIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let lz = self.set.leading_zeros();
        if lz == self.set.len() {
            None
        } else {
            self.set.remove_one(lz);
            Some(lz)
        }
    }
}

impl IndexMut<usize> for Graph {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl Index<usize> for Graph {
    type Output = Set;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

fn bitvec_from_f(n: usize, i: usize, f: fn(usize, usize, usize) -> bool) -> BitVec<usize, Msb0> {
    let mut result: BitVec<usize, Msb0> = bitvec![usize, Msb0; 0; n];
    for i_other_vertex in 0..n {
        if i_other_vertex == i {
            assert!(
                !f(i, i_other_vertex, n),
                "function given to from_f must not be true for i == i_other_vertex"
            )
        }
        result.set(i_other_vertex, f(i, i_other_vertex, n));
    }
    result
}

fn bitvec_from_closure<F: Fn(usize, usize, usize) -> bool + Copy>(
    n: usize,
    i: usize,
    f: F,
) -> BitVec<usize, Msb0> {
    let mut result: BitVec<usize, Msb0> = bitvec![usize, Msb0; 0; n];
    for i_other_vertex in 0..n {
        if i_other_vertex == i {
            assert!(
                !f(i, i_other_vertex, n),
                "closure given to from_closure must not be true for i == i_other_vertex"
            )
        }
        result.set(i_other_vertex, f(i, i_other_vertex, n));
    }
    result
}

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
        todo!()
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

pub fn u32_to_bitvec(u: u32, n: usize) -> BitVec<usize, Msb0> {
    let mut bv: BitVec<usize, Msb0> = bitvec![usize, Msb0;];
    bv.extend_from_bitslice(u.view_bits::<Msb0>().split_at(n).0);
    bv
}

/// a dynamic associative array usize -> usize based on Vec rather than HashMap or BTreeMap
struct VecMap(Vec<Option<usize>>);

impl VecMap {
    pub fn new() -> Self {
        Self(vec![])
    }

    pub fn get(&self, index: usize) -> usize {
        self.get_safely(index).unwrap()
    }

    pub fn get_safely(&self, index: usize) -> Option<usize> {
        if index >= self.0.len() {
            None
        } else {
            self.0[index]
        }
    }

    pub fn set(&mut self, index: usize, value: usize) {
        while self.0.len() < index + 1 {
            self.0.push(None);
        }
        self.0[index] = Some(value);
    }

    pub fn increment(&mut self, index: usize) {
        if let Some(v) = &mut self.0[index] {
            *v += 1;
        }
    }
}

// static variables in nauty.c
#[derive(Default)]
pub struct NautyEnv {
    /* temporary versions of some stats: */
    pub invapplics: usize,
    pub invsuccesses: usize,
    pub invarsuclevel: usize,
    pub noncheaplevel: usize, /* level of greatest ancestor for which cheapautom==FALSE */
    pub eqlev_canon: isize,   /* level to which codes for this node match those for the bsf leaf. */

    pub needshortprune: bool, /* used to flag calls to shortprune */

    pub workperm: Vec<usize>,
    pub active: Vec<Set>,
    pub workspace: Vec<Set>, /*work area to hold automorphism data */
}

/*****************************************************************************
*                                                                            *
*  This procedure finds generators for the automorphism group of a           *
*  vertex-coloured graph and optionally finds a canonically labelled         *
*  isomorph.  A description of the data structures can be found in           *
*  nauty.h and in the "nauty User's Guide".  The Guide also gives            *
*  many more details about its use, and implementation notes.                *
*                                                                            *
*  Parameters - <r> means read-only, <w> means write-only, <wr> means both:  *
*           g <r>  - the graph                                               *
*     lab,ptn <rw> - used for the partition nest which defines the colouring *
*                  of g.  The initial colouring will be set by the program,  *
*                  using the same colour for every vertex, if                *
*                  options->defaultptn!=FALSE.  Otherwise, you must set it   *
*                  yourself (see the Guide). If options->getcanon!=FALSE,    *
*                  the contents of lab on return give the labelling of g     *
*                  corresponding to canong.  This does not change the        *
*                  initial colouring of g as defined by (lab,ptn), since     *
*                  the labelling is consistent with the colouring.           *
*     active  <r>  - If this is not NULL and options->defaultptn==FALSE,     *
*                  it is a set indicating the initial set of active colours. *
*                  See the Guide for details.                                *
*     orbits  <w>  - On return, orbits[i] contains the number of the         *
*                  least-numbered vertex in the same orbit as i, for         *
*                  i=0,1,...,n-1.                                            *
*    options  <r>  - A list of options.  See nauty.h and/or the Guide        *
*                  for details.                                              *
*      stats  <w>  - A list of statistics produced by the procedure.  See    *
*                  nauty.h and/or the Guide for details.                     *
*  workspace  <w>  - A chunk of memory for working storage.                  *
*  worksize   <r>  - The number of setwords in workspace.  See the Guide     *
*                  for guidance.                                             *
*          m  <r>  - The number of setwords in sets.  This must be at        *
*                  least ceil(n / WORDSIZE) and at most MAXM.                *
*          n  <r>  - The number of vertices.  This must be at least 1 and    *
*                  at most MAXN.                                             *
*     canong  <w>  - The canononically labelled isomorph of g.  This is      *
*                  only produced if options->getcanon!=FALSE, and can be     *
*                  given as NULL otherwise.                                  *
*                                                                            *
*  FUNCTIONS CALLED: firstpathnode(),updatecan()                             *
*                                                                            *
*****************************************************************************/
fn nauty(
    g_arg: Graph,
    lab: &mut [usize],
    ptn: &mut [usize],
    active_arg: &[Set],
    orbits_arg: &mut Vec<usize>,
    options: &OptionBlk,
    stats_arg: &mut StatBlk,
    canong_arg: &mut Graph,
) -> Result<(), u8> {
    let mut nauty_env = NautyEnv::default();
    let n = g_arg.n();

    let mut numcells: usize;
    let mut initstatus: u8;

    let defltwork: Vec<Set>;
    let workperm: Vec<usize>;
    let fixedpts: Vec<Set>;
    let firstlab: Vec<usize>;
    let canonlab: Vec<usize>;
    let firstcode: Vec<u8>;
    let canoncode: Vec<u8>;
    let firsttc: Vec<usize>;
    let mut active: Vec<Set>;

    /* initialize everything: */
    if options.defaultptn {
        for i in 0..n {
            lab[i] = i;
            ptn[i] = NAUTY_INFINITY;
        }
        ptn[n - 1] = 0;
        active = vec![(Set::new())];
        active[0].add_one(0);
        numcells = 1;
    } else {
        ptn[n - 1] = 0;
        numcells = 0;
        for i in 0..n {
            if ptn[i] != 0 {
                ptn[i] = NAUTY_INFINITY;
            } else {
                numcells += 1;
            }
            if active_arg.is_empty() {
                active = vec![(Set::new())];
                for mut i in 0..n {
                    active[0].add_one(i);
                    while ptn[i] != 0 {
                        i += 1;
                    }
                }
            } else {
                active = active_arg.to_vec();
            }
        }
    }
    let mut g: Graph;
    let mut cannong: Graph;
    initstatus = 0;

    let mut orbits: Vec<usize> = (0..n).collect();
    let mut stats: StatBlk = StatBlk::new(n);
    fixedpts = vec![];
    nauty_env.noncheaplevel = 1;
    nauty_env.eqlev_canon = -1;
    nauty_env.needshortprune = false;
    nauty_env.invarsuclevel = NAUTY_INFINITY;
    nauty_env.invapplics = 0;
    nauty_env.invsuccesses = 0;
    firstpathnode0(g_arg, lab, ptn, 1, numcells, LinkedList::new(), &mut stats);
    Ok(())
}

/*****************************************************************************
*                                                                            *
*  firstpathnode(lab,ptn,level,numcells) produces a node on the leftmost     *
*  path down the tree.  The parameters describe the level and the current    *
*  colour partition.  The set of active cells is taken from the global set   *
*  'active'.  If the refined partition is not discrete, the leftmost child   *
*  is produced by calling firstpathnode, and the other children by calling   *
*  othernode.                                                                *
*  For MAXN=0 there is an extra parameter: the address of the parent tcell   *
*  structure.                                                                *
*  The value returned is the level to return to.                             *
*                                                                            *
*  FUNCTIONS CALLED: (*usernodeproc)(),doref(),cheapautom(),                 *
*                    firstterminal(),nextelement(),breakout(),               *
*                    firstpathnode(),othernode(),recover(),writestats(),     *
*                    (*userlevelproc)(),(*tcellproc)(),shortprune()          *
*                                                                            *
*****************************************************************************/
fn firstpathnode0(
    g_arg: Graph,
    lab: &mut [usize],
    ptn: &mut [usize],
    level: usize,
    numcells: usize,
    tcnode_parent: LinkedList<Set>,
    stats: &mut StatBlk,
) -> Result<(), u8> {
    let tv: usize;
    let tv1: usize;
    let index: usize;
    let rtnlevel: usize;
    let tcellsize: usize;
    let tc: usize;
    let childcount: usize;
    let qinvar: usize;
    let refcode: usize;
    let mut tcell: &mut Set;
    let mut tcnode_this: LinkedList<Set> = tcnode_parent;

    if tcnode_this.is_empty() {
        tcnode_this.push_back(Set::new());
    }
    tcell = tcnode_this.front_mut().unwrap();
    stats.numnodes += 1;

    /* refine partition : */
    // doref(
    //     g_arg,
    //     lab,
    //     ptn,
    //     level,
    //     &numcells,
    //     &qinvar,
    //     workperm,
    //     active,
    //     &refcode,
    //     dispatch.refine,
    //     invarproc,
    //     mininvarlevel,
    //     maxinvarlevel,
    //     invararg,
    //     digraph,
    //     M,
    //     n,
    // );
    Ok(())
}
