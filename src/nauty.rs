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
    naugraph::isautom,
    nautil::{doref_nest, maketargetcell, maketargetcell_mut},
    nauty::partition_nest::{
        PartitionNest,
        partition::{Partition, cell::Cell},
    },
};
use bitvec::{bitvec, order::Msb0, vec::BitVec, view::BitView};
use std::{
    fmt::Debug,
    fs::File,
    mem,
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
    tc_level: usize,    /* max level for smart target cell choosing */
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
pub mod graph;
pub type NautyCounter = u128;

pub trait SetTrait {
    fn difference(&self, other: &Self) -> Self;
    fn first_bit_nz_index(&self) -> usize;
    fn zeros(n: usize) -> Self;
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

    fn zeros(n: usize) -> Self {
        bitvec![usize, Msb0; 0; n]
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

impl<'a, 'b> From<&'a Cell<'b>> for Set {
    fn from(value: &'a Cell) -> Self {
        let mut set = Self::zeros(value.partition.len());
        for i in value.iter() {
            set.set(*i, true);
        }
        set
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

impl IndexMut<usize> for graph::Graph {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl Index<usize> for graph::Graph {
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
    /* working variables: <the "bsf leaf" is the leaf which is best guess so far at the canonical leaf>  */
    pub gca_first: usize, /* level of greatest common ancestor of current node and first leaf */
    pub gca_canon: usize, /* ditto for current node and bsf leaf */
    pub noncheaplevel: usize, /* level of greatest ancestor for which cheapautom==FALSE */
    pub allsamelevel: usize, /* level of least ancestor of first leaf for which all descendant leaves are known to be equivalent */
    pub eqlev_first: usize,  /* level to which codes for this node match those for first leaf */
    pub eqlev_canon: isize,  /* level to which codes for this node match those for the bsf leaf. */
    pub comp_canon: isize, /* -1,0,1 according as code at eqlev_canon+1 is <,==,> that for bsf leaf.  Also used for similar purpose during leaf processing */
    pub samerows: usize, /* number of rows of canong which are correct for the bsf leaf  BDM:correct description? */
    pub canonlevel: usize, /* level of bsf leaf */
    pub stabvertex: usize, /* point fixed in ancestor of first leaf at level gca_canon */
    pub cosetindex: usize, /* the point being fixed at level gca_first */
    pub needshortprune: bool, /* used to flag calls to shortprune */

    pub workperm: Vec<usize>,
    pub fixedpts: Set,
    pub first_partition: Partition,
    pub canon_partition: Partition,
    pub firstcode: Vec<u16>,
    pub canoncode: Vec<u16>,
    pub firsttc: Vec<isize>,
    pub active: Vec<Set>,
    pub workspace: Vec<Set>, /*work area to hold automorphism data */
}

impl NautyEnv {
    fn new(n: usize) -> Self {
        Self {
            firstcode: vec![0; n + 2],
            firsttc: vec![0; n + 2],
            ..Default::default()
        }
    }
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
    g_arg: graph::Graph,
    lab: &mut [usize],
    ptn: &mut [usize],
    active_arg: &Set,
    orbits_arg: &mut Vec<usize>,
    options: &OptionBlk,
    _stats_arg: &mut StatBlk,
    _canong_arg: &mut graph::Graph,
) -> Result<(), u8> {
    let n = g_arg.n();
    let mut nauty_env = NautyEnv::new(n);

    let mut numcells: usize;
    let mut _initstatus: u8;

    let _defltwork: Vec<Set>;
    let mut firstcode: Vec<usize> = vec![0; n + 2];
    let _canoncode: Vec<u8>;
    let mut firsttc: VecMap = VecMap::new();
    let mut active: Set;

    /* initialize everything: */
    active = Set::one(0);
    if options.defaultptn {
        for i in 0..n {
            lab[i] = i;
            ptn[i] = NAUTY_INFINITY;
        }
        ptn[n - 1] = 0;
        active = Set::one(0);
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
                active = Set::zeros(n);
                for mut i in 0..n {
                    active.add_one(i);
                    while ptn[i] != 0 {
                        i += 1;
                    }
                }
            } else {
                active = active_arg.clone();
            }
        }
    }
    let mut g: graph::Graph;
    let mut cannong: graph::Graph;
    _initstatus = 0;

    *orbits_arg = (0..n).collect();
    let mut stats: StatBlk = StatBlk::new(n);
    nauty_env.fixedpts = Set::zeros(n);
    nauty_env.noncheaplevel = 1;
    nauty_env.eqlev_canon = -1;
    nauty_env.needshortprune = false;
    nauty_env.invarsuclevel = NAUTY_INFINITY;
    nauty_env.invapplics = 0;
    nauty_env.invsuccesses = 0;
    let mut partition = Partition::new(PartitionNest::new(lab.to_vec(), ptn.to_vec()), 1);
    firstpathnode_nest(
        &g_arg,
        &mut partition,
        &mut active,
        &mut firstcode,
        &mut stats,
        options.tc_level,
        &mut firsttc,
        orbits_arg,
        &mut nauty_env,
        &options,
    );
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
fn firstpathnode(
    _g_arg: graph::Graph,
    _lab: &mut [usize],
    _ptn: &mut [usize],
    _level: usize,
    _numcells: usize,
    stats: &mut StatBlk,
    _firsttc: &mut VecMap,
) -> Result<(), u8> {
    let _tv: usize;
    let _tv1: usize;
    let _index: usize;
    let _rtnlevel: usize;
    let _tcellsize: usize;
    let _tc: usize;
    let _childcount: usize;
    let _qinvar: usize;
    let _refcode: usize;
    let mut _tcell: &mut Set;

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
// 561
fn firstpathnode_nest(
    mut g_arg: &graph::Graph,
    partition: &mut Partition,
    mut active: &mut Set,
    firstcode: &mut Vec<usize>,
    stats: &mut StatBlk,
    tc_level: usize,
    firsttc: &mut VecMap,
    orbits_arg: &mut Vec<usize>,
    nauty_env: &mut NautyEnv,
    options: &OptionBlk,
) -> usize {
    let _index: usize;
    let mut rtnlevel: usize = 0;
    let _childcount: usize;
    let mut qinvar: usize = 0;
    let mut refcode: usize = 0;
    let mut partition = partition;
    let level = partition.level;
    let cheapautom = partition.cheapautom();

    stats.numnodes += 1;

    /* refine partition : */
    doref_nest(
        g_arg,
        &mut partition,
        &mut qinvar,
        &mut active,
        &mut refcode,
    );
    firstcode[partition.level] = refcode;
    if qinvar > 0 {
        todo!("qinvar always == 0");
    }
    if !partition.is_discrete() {
        let tcell = maketargetcell(&g_arg, &mut partition, tc_level, None);
        stats.tctotal += tcell.len();
        nauty_env.firsttc[level] = tcell.first_lab_index as isize;
    } else {
        nauty_env.firsttc[level] = -1;
    }
    if partition.is_discrete() {
        firstterminal(&mut partition, stats, nauty_env, options.getcanon);
        return partition.level - 1;
    }

    let tcell = maketargetcell(&g_arg, &mut partition, tc_level, None);
    stats.tctotal += tcell.len();
    firsttc.set(level, tcell.first_lab_index);
    if nauty_env.noncheaplevel >= level && !cheapautom {
        nauty_env.noncheaplevel += 1;
    }

    /* use the elements of the target cell to produce the children: */
    let mut _index = 0;
    let mut _childcount = 0;
    let cell_lab = tcell.to_vec();
    let tv1 = tcell[0];
    let tc = tcell.first_lab_index;
    for (i, tv) in cell_lab.iter().copied().enumerate() {
        if orbits_arg[tv] == tv {
            partition.breakout(tc, tv);
            nauty_env.fixedpts.add_one(tv);
            nauty_env.cosetindex = tv;
            if tv == tv1 {
                partition.advance();
                rtnlevel = firstpathnode_nest(
                    g_arg, partition, active, firstcode, stats, tc_level, firsttc, orbits_arg,
                    nauty_env, options,
                );
                _childcount = 1;
                nauty_env.gca_first = level;
                nauty_env.stabvertex = tv1;
            } else {
                partition.advance();
                //                rtnlevel = othernode(partition);
                _childcount += 1;
            }
            if rtnlevel < level {
                return rtnlevel;
            }
            if nauty_env.needshortprune {
                //                shortprune(cell_lab)
            }
            //            recover(partition.nest, level)
        }
    }

    0
}

/*****************************************************************************
*                                                                            *
*  othernode(lab,ptn,level,numcells) produces a node other than an ancestor  *
*  of the first leaf.  The parameters describe the level and the colour      *
*  partition.  The list of active cells is found in the global set 'active'. *
*  The value returned is the level to return to.                             *
*                                                                            *
*  FUNCTIONS CALLED: (*usernodeproc)(),doref(),refine(),recover(),           *
*                    processnode(),cheapautom(),(*tcellproc)(),shortprune(), *
*                    nextelement(),breakout(),othernode(),longprune()        *
*                                                                            *
*****************************************************************************/
// 720
fn othernode(
    mut g_arg: &graph::Graph,
    partition: &mut Partition,
    mut active: &mut Set,
    stats: &mut StatBlk,
    tc_level: usize,
    firsttc: &mut VecMap,
    orbits_arg: &mut Vec<usize>,
    nauty_env: &mut NautyEnv,
    options: &OptionBlk,
) -> usize {
    let _index: usize;
    let mut _rtnlevel: usize;
    let _childcount: usize;
    let mut qinvar: usize = 0;
    let mut refcode: usize = 0;
    let mut code: u16 = 0;
    let mut partition = partition;
    let level = partition.level;
    let _cheapautom = partition.cheapautom();

    stats.numnodes += 1;
    /* refine partition : */
    doref_nest(
        g_arg,
        &mut partition,
        &mut qinvar,
        &mut active,
        &mut refcode,
    );
    code = refcode as u16;

    if qinvar > 0 {
        nauty_env.invapplics += 1;
        if qinvar == 2 {
            nauty_env.invsuccesses += 1;
            if level < nauty_env.invarsuclevel {
                nauty_env.invarsuclevel = level;
            }
        }
    }

    if nauty_env.eqlev_first == level - 1 && code == nauty_env.firstcode[level] {
        nauty_env.eqlev_first = level;
    }
    if options.getcanon != 0 {
        if nauty_env.eqlev_canon == level as isize - 1 {
            if code < nauty_env.canoncode[level] {
                nauty_env.comp_canon = -1;
            } else if code > nauty_env.canoncode[level] {
                nauty_env.comp_canon = 1;
            } else {
                nauty_env.comp_canon = 0;
                nauty_env.eqlev_canon = level as isize;
            }
        }
        if nauty_env.comp_canon > 0 {
            nauty_env.canoncode[level] = code;
        }
    }

    let tc: isize = -1;
    if !partition.is_discrete()
        && (nauty_env.eqlev_first == level || (options.getcanon != 0 && nauty_env.comp_canon >= 0))
    {
        let tcell;
        if options.getcanon == 0 || nauty_env.comp_canon < 0 {
            tcell = maketargetcell(&g_arg, &mut partition, tc_level, Some(firsttc.get(level)));
            if tc != firsttc.get(level) as isize {
                nauty_env.eqlev_first = level - 1;
            }
        } else {
            tcell = maketargetcell(&g_arg, &mut partition, tc_level, None);
        }
        stats.tctotal += tcell.len();
    }

    0
}

/*****************************************************************************
*                                                                            *
*  Process the first leaf of the tree.                                       *
*                                                                            *
*  FUNCTIONS CALLED: NONE                                                    *
*                                                                            *
*****************************************************************************/
//866

fn firstterminal(
    partition: &Partition,
    stats: &mut StatBlk,
    nauty_env: &mut NautyEnv,
    getcanon: u8,
) {
    let level = partition.level;
    stats.maxlevel = level;
    nauty_env.gca_first = level;
    nauty_env.allsamelevel = level;
    nauty_env.eqlev_first = level;
    nauty_env.firstcode[level + 1] = 0o77777;
    nauty_env.firsttc[level + 1] = -1;

    nauty_env.first_partition = partition.clone();

    if getcanon != 0 {
        nauty_env.canonlevel = level;
        nauty_env.eqlev_canon = level as isize;
        nauty_env.gca_canon = level;
        nauty_env.comp_canon = 0;
        nauty_env.samerows = 0;
        nauty_env.canon_partition = partition.clone();
        nauty_env.canoncode[0..=level].copy_from_slice(&nauty_env.firstcode[0..=level]);
        nauty_env.canoncode[level + 1] = 0o77777;
        stats.canupdates = 1;
    }
}

/*****************************************************************************
*                                                                            *
*  Process a node other than the first leaf or its ancestors.  It is first   *
*  classified into one of five types and then action is taken appropriate    *
*  to that type.  The types are                                              *
*                                                                            *
*  0:   Nothing unusual.  This is just a node internal to the tree whose     *
*         children need to be generated sometime.                            *
*  1:   This is a leaf equivalent to the first leaf.  The mapping from       *
*         firstlab to lab is thus an automorphism.  After processing the     *
*         automorphism, we can return all the way to the closest invocation  *
*         of firstpathnode.                                                  *
*  2:   This is a leaf equivalent to the bsf leaf.  Again, we have found an  *
*         automorphism, but it may or may not be as useful as one from a     *
*         type-1 node.  Return as far up the tree as possible.               *
*  3:   This is a new bsf node, provably better than the previous bsf node.  *
*         After updating canonlab etc., treat it the same as type 4.         *
*  4:   This is a leaf for which we can prove that no descendant is          *
*         equivalent to the first or bsf leaf or better than the bsf leaf.   *
*         Return up the tree as far as possible, but this may only be by     *
*         one level.                                                         *
*                                                                            *
*  Types 2 and 3 can't occur if getcanon==FALSE.                             *
*  The value returned is the level in the tree to return to, which can be    *
*  anywhere up to the closest invocation of firstpathnode.                   *
*                                                                            *
*  FUNCTIONS CALLED:    isautom(),updatecan(),testcanlab(),fmperm(),         *
*                       writeperm(),(*userautomproc)(),orbjoin(),            *
*                       shortprune(),fmptn()                                 *
*                                                                            *
*****************************************************************************/
// 923
fn processnode(
    g: &graph::Graph,
    partition: &Partition,
    nauty_env: &mut NautyEnv,
    options: &OptionBlk,
) -> usize {
    let mut code: u8 = 0;
    let mut newlevel: usize = 0;

    if nauty_env.eqlev_first != partition.level
        && (options.getcanon == 0 || nauty_env.comp_canon < 0)
    {
        code = 4;
    } else if partition.is_discrete() {
        if nauty_env.eqlev_first == partition.level {
            let perm = partition.permutation(&nauty_env.first_partition);
            if nauty_env.gca_first >= nauty_env.noncheaplevel || isautom(g, &perm) {
                code = 1;
            }
        }
    }
    if code == 0 {
        if options.getcanon != 0 {
            let mut sr = 0;
            if nauty_env.comp_canon == 0 {
                if partition.level < nauty_env.canonlevel {
                    nauty_env.comp_canon = 1;
                } else {
                }
            }
        }
    }

    newlevel
}
