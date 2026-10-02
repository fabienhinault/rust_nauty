use crate::{
    naugraph::{refine_nest, targetcell, targetcell_mut},
    nauty::{
        NautyEnv, Set, SetTrait,
        graph::{self, Graph},
        partition_nest::partition::{Partition, cell::Cell, cell_mut::CellMut},
    },
};
use std::{
    fs::File,
    io::{BufWriter, Write},
};

pub trait SetWordNautilTrait {
    fn next_element(&self, pos: Option<usize>) -> Option<usize>;
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
    // 111
    // should be generally replaced by iteration on Set
    fn next_element(&self, pos: Option<usize>) -> Option<usize> {
        let setwd: Set = match pos {
            Some(pos) => self.filter(&self.bit_mask(pos)),
            None => self.clone(),
        };
        let lz = setwd.leading_zeros();
        if lz == setwd.len() { None } else { Some(lz) }
    }
}

/*****************************************************************************
*                                                                            *
*  orbits represents a partition of {0,1,...,n-1}, by orbits[i] = the        *
*  smallest element in the same cell as i.  map[] is any array with values   *
*  in {0,1,...,n-1}.  orbjoin(orbits,map,n) joins the cells of orbits[]      *
*  together to the minimum extent such that for each i, i and map[i] are in  *
*  the same cell.  The function value returned is the new number of cells.   *
*                                                                            *
*  GLOBALS ACCESSED: NONE                                                    *
*                                                                            *
*****************************************************************************/
// 266
pub fn orbjoin(orbits: &mut [usize], map: &[usize]) -> usize {
    let mut j1;
    let mut j2;
    for i in 0..orbits.len() {
        if map[i] != i {
            j1 = orbits[i];
            while orbits[j1] != j1 {
                j1 = orbits[j1];
            }
            j2 = orbits[map[i]];
            while orbits[j2] != j2 {
                j2 = orbits[j2];
            }
            if j1 < j2 {
                orbits[j2] = j1;
            } else if j1 > j2 {
                orbits[j1] = j2;
            }
        }
    }
    for i in 0..orbits.len() {
        orbits[i] = orbits[orbits[i]];
    }
    orbits
        .iter()
        .enumerate()
        .filter(|(i, orbit)| *orbit == i)
        .count()
}

fn condnl(
    x: usize,
    curlen: &mut usize,
    linelength: Option<usize>,
    writer: &mut BufWriter<&mut File>,
) {
    if let Some(linelength) = linelength
        && *curlen + x > linelength
    {
        write!(writer, "\n   ");
        *curlen = 3;
    }
}

/*****************************************************************************
*                                                                            *
*  writeperm(f,perm,cartesian,linelength,n) writes the permutation perm to   *
*  the file f.  The cartesian representation (i.e. perm itself) is used if   *
*  cartesian != FALSE; otherwise the cyclic representation is used.  No      *
*  more than linelength characters (not counting '\n') are written on each   *
*  line, unless linelength is ridiculously small.  linelength<=0 causes no   *
*  line breaks at all to be made.  The global int labelorg is added to each  *
*  vertex number.                                                            *
*                                                                            *
*  GLOBALS ACCESSED: itos(),putstring()                                      *
*                                                                            *
*****************************************************************************/
// 304
pub fn writeperm(
    f: &mut File,
    perm: &[usize],
    cartesian: bool,
    linelength: Option<usize>,
    nauty_env: &NautyEnv,
) {
    let mut writer = BufWriter::new(f);
    let labelorg = nauty_env.labelorg;
    if cartesian {
        let mut curlen = 0;
        for p in perm {
            let s = format!("{}", p + labelorg);
            condnl(s.len() + 1, &mut curlen, linelength, &mut writer);
            write!(writer, " ").expect("write!");
            write!(writer, "{s}").expect("write!");
            curlen += s.len() + 1;
        }
        writeln!(writer).expect("writeln!");
    } else {
        let mut curlen = 0;
        let mut workperm = vec![false; perm.len()];
        for i in 0..perm.len() {
            if !workperm[i] && perm[i] != i {
                let mut l = i;
                let mut s = format!("{}", l + labelorg);
                if curlen > 3 {
                    condnl(2 * s.len() + 4, &mut curlen, linelength, &mut writer);
                    write!(writer, "(").expect("write!");
                    loop {
                        write!(writer, "{s}").expect("write!");
                        curlen += s.len() + 1;
                        let k = l;
                        l = perm[l];
                        workperm[k] = true;
                        if l != i {
                            s = format!("{}", l + labelorg);
                            condnl(s.len(), &mut curlen, linelength, &mut writer);
                            write!(writer, " ").expect("write!");
                        }
                        if l == i {
                            break;
                        }
                    }
                    write!(writer, ")").expect("write!");
                    curlen += 1;
                }
            }
        }
        if curlen == 0 {
            writeln!(writer, "(1)").expect("write!");
        } else {
            writeln!(writer).expect("write!");
        }
    }
}

/*****************************************************************************
*                                                                            *
*  fmperm(perm,fix,mcr,m,n) uses perm to construct fix and mcr.  fix         *
*  contains those points are fixed by perm, while mcr contains the set of    *
*  those points which are least in their orbits.                             *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>                                                  *
*                                                                            *
*****************************************************************************/
// 379
pub fn fmperm(perm: &[usize]) -> (Set, Set) {
    let mut fix = Set::zeros(perm.len());
    let mut mcr = Set::zeros(perm.len());
    let mut workperm = vec![0; perm.len()];
    for i in 0..perm.len() {
        if perm[i] == i {
            fix.set(i, true);
            mcr.set(i, true);
        } else if workperm[i] == 0 {
            let mut l = i;
            let mut k;
            loop {
                k = l;
                l = perm[l];
                workperm[k] = 1;
                if l != i {
                    break;
                }
            }
            mcr.set(i, true);
        }
    }
    (fix, mcr)
}

/*****************************************************************************
*                                                                            *
*  doref(g,lab,ptn,level,numcells,qinvar,invar,active,code,refproc,          *
*        invarproc,mininvarlev,maxinvarlev,invararg,digraph,m,n)             *
*  is used to perform a refinement on the partition at the given level in    *
*  (lab,ptn).  The number of cells is *numcells both for input and output.   *
*  The input active is the active set for input to the refinement procedure  *
*  (*refproc)(), which must have the argument list of refine().              *
*  active may be arbitrarily changed.  invar is used for working storage.    *
*  First, (*refproc)() is called.  Then, if invarproc!=NULL and              *
*  |mininvarlev| <= level <= |maxinvarlev|, the routine (*invarproc)() is    *
*  used to compute a vertex-invariant which may refine the partition         *
*  further.  If it does, (*refproc)() is called again, using an active set   *
*  containing all but the first fragment of each old cell.  Unless g is a    *
*  digraph, this guarantees that the final partition is equitable.  The      *
*  arguments invararg and digraph are passed to (*invarproc)()               *
*  uninterpretted.  The output argument code is a composite of the codes     *
*  from all the calls to (*refproc)().  The output argument qinvar is set    *
*  to 0 if (*invarproc)() is not applied, 1 if it is applied but fails to    *
*  refine the partition, and 2 if it succeeds.                               *
*  See the file nautinv.c for a further discussion of vertex-invariants.     *
*  Note that the dreadnaut I command generates a call to  this procedure     *
*  with level = mininvarlevel = maxinvarlevel = 0.                           *
*                                                                            *
*****************************************************************************/
// 481
// case where invarproc is null, dorest just calls refine
pub fn doref_nest(
    g: &graph::Graph,
    partition: &mut Partition,
    qinvar: &mut usize,
    active: &mut Set,
    code: &mut usize,
) {
    refine_nest(g, partition, active, code);
    *qinvar = 0;
}

/*****************************************************************************
*                                                                            *
*  maketargetcell(g,lab,ptn,level,tcell,tcellsize,&cellpos,                  *
*                 tc_level,digraph,hint,targetcell,m,n)                      *
*  calls targetcell() to determine the target cell at the specified level    *
*  in the partition nest (lab,ptn).  It must be a nontrivial cell (if not,   *
*  the first cell.  The intention of hint is that, if hint >= 0 and there    *
*  is a suitable non-trivial cell starting at position hint in lab,          *
*  that cell is chosen.                                                      *
*  tc_level and digraph are input options.                                   *
*  When a cell is chosen, tcell is set to its contents, *tcellsize to its    *
*  size, and cellpos to its starting position in lab.                        *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>                                                  *
*                                                                            *
*****************************************************************************/
// l565
pub fn maketargetcell<'a>(
    g: &graph::Graph,
    partition: &'a Partition,
    tc_level: usize,
    hint: Option<usize>,
) -> Cell<'a> {
    targetcell(g, partition, hint, tc_level)
}

pub fn maketargetcell_mut<'a>(
    g: &graph::Graph,
    partition: &'a mut Partition,
    tc_level: usize,
    hint: Option<usize>,
) -> CellMut<'a> {
    targetcell_mut(g, partition, hint, tc_level)
}
