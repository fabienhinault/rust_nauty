use crate::{
    naugraph::refine_nest,
    nauty::{Graph, Set, SetTrait, partition_nest::partition::Partition},
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

// pub fn doref(
//     g: Graph,
//     lab: &mut [usize],
//     ptn: &mut [usize],
//     level: usize,
//     numcells: &mut usize,
//     qinvar: &mut usize,
//     invar: &mut Vec<usize>,
//     active: &mut Set,
//     code: &mut usize,
//     refproc: RP,
//     invarproc: IP,
//     mininvarlev: usize,
//     maxinvarlev: usize,
//     invararg: usize,
//     digraph: bool,
//     nauty_env: &mut NautyEnv,
// ) {
//     let pw: usize;
//     let i: usize;
//     let cell1: usize;
//     let cell2: usize;
//     let nc: usize;
//     let tvpos: usize;
//     let minlev: usize;
//     let maxlev: usize;
//     let longcode: usize;
//     let same: bool;
//     nauty_env.workperm = Vec::with_capacity(g.n());

//     tvpos = active.first_one().unwrap_or(0);

//     refproc()
// }

// case where invarproc is null, dorest just calls refine
pub fn doref_nest(
    g: &mut Graph,
    partition: &mut Partition,
    qinvar: &mut usize,
    active: &mut Set,
    code: &mut usize,
) {
    refine_nest(g, partition, active, code);
    *qinvar = 0;
}

pub fn maketargetcell(
    g: &Graph,
    partition: &Partition,
    tcell: &mut Set,
    tc_level: usize,
    hint: usize,
) {
    let i: usize;
    let j: usize;
    let k: usize;
}
