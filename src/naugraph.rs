use crate::{
    nautil::SetWordNautilTrait,
    nauty::{Graph, NautyEnv, Set, SetTrait},
};
use cfor::cfor;
use core::num;

struct NaugraphEnv {
    pub workset: Vec<Set>,
    pub workperm: Vec<usize>,
    pub bucket: Vec<usize>,
    pub dnwork: Vec<Set>,
}

/* macros for hash-codes: */
/* : expression whose long value depends only on long l and int/long i.
Anything goes, preferably non-commutative. */
fn mash(l: usize, i: usize) -> usize {
    ((l ^ 0o65435) + i) & 0o77777
}
fn cleanup(l: usize) -> usize {
    l % 0o77777
}

/*****************************************************************************
*                                                                            *
*  isautom(g,perm,digraph,m,n) = TRUE iff perm is an automorphism of g       *
*  (i.e., g^perm = g).  Symmetry is assumed unless digraph = TRUE.           *
*                                                                            *
*****************************************************************************/
fn is_autom(g: &Graph, perm: &[usize]) -> bool {
    for (i_row, row) in g.0.iter().enumerate() {
        let p_row = g.get(perm[i_row]);
        for pos in row.masked(i_row).ones_iter() {
            if !p_row[perm[pos]] {
                return false;
            }
        }
    }
    true
}
///
///
///    g: &mut Graph,   
///    lab: &mut [usize],        labels
///    ptn: &mut [isize],        partition
///    level: isize,             recursion level
///    numcells: &mut usize,     number of cells
///    count: &mut Vec<usize>,   number of vertices in cells
///    active: &mut Set,         
///    code: &mut usize,         
fn refine(
    g: &mut Graph,
    lab: &mut [usize],
    ptn: &mut [isize],
    level: isize,
    numcells: &mut usize,
    count: &mut Vec<usize>,
    active: &mut Set,
    code: &mut usize,
) {
    let mut i: usize;
    let mut c1: usize;
    let mut c2: usize;
    let mut labc1: usize;
    let mut split1: usize;
    let mut split2: usize;
    let mut cell1: usize;
    let mut cell2: usize;
    let mut cnt: usize;
    let mut bmin: usize;
    let mut bmax: usize;
    let mut longcode: usize;
    let mut gptr: &Set;
    let mut maxcell: isize;
    let mut maxpos: Option<usize> = None;
    let mut hint: usize;

    let workperm: Vec<usize> = Vec::with_capacity(g.n());
    let mut bucket: Vec<usize> = Vec::with_capacity(g.n() + 2);

    longcode = *numcells;
    hint = 0;
    loop {
        if *numcells == g.n() {
            break;
        }
        split1 = hint;
        if !active[split1] {
            split1 = match active.next_element(Some(split1)) {
                Some(next) => next,
                None => match active.first_one() {
                    Some(first) => first,
                    None => break,
                },
            }
        }
        active.remove_one(split1);
        split2 = split1;
        while ptn[split2] > level {
            split2 += 1;
        }
        longcode = mash(longcode, split1 + split2);
        /* trivial splitting cell */
        if split1 == split2 {
            gptr = &g.0[lab[split1]];
            cell2 = 0;
            cfor! {cell1 = 0; cell1 < g.n(); cell1 = cell2 + 1; {
                cfor! {cell2 = cell1; ptn[cell2] > level; cell2 += 1; {}}
                if cell1 == cell2 {
                    continue;
                }
                c1 = cell1;
                c2 = cell2;
                while c1 <= c2 {
                    labc1 = lab[c1];
                    if gptr[labc1] {
                        c1 +=1;
                    } else {
                        lab[c1] = lab[c2];
                        lab[c2] = labc1;
                        c2 -= 1;
                    }
                }
                if c2 >= cell1 && c1 <= cell2{
                    ptn[c2] = level;
                    longcode = mash(longcode, c2);
                    *numcells += 1;
                    if active[cell1] || (c2-cell1 >= cell2-c1) {
                        active.add_one(c1);
                        if c1 == cell2 {
                            hint = c1;
                        }
                    } else {
                        active.add_one(cell1);
                        if c2 == cell1 {
                            hint = cell1;
                        }
                    }
                }
            }}
        /* nontrivial splitting cell */
        } else {
            let mut workset: Set = Set::new();
            for i in split1..=split2 {
                workset.add_one(lab[i]);
            }
            longcode = mash(longcode, split2 - split1 + 1);

            cell2 = 0;
            cfor! {cell1 = 0; cell1 < g.n(); cell1 = cell2 + 1;
            {
                cfor! {cell2 = cell1; ptn[cell2] > level; cell2 += 1; {}}
                if cell1 == cell2 {
                    continue;
                }
                i = cell1;
                cnt = (workset.clone() & g.0[lab[i]].clone()).count_ones();
                count[i] = cnt;
                bmin = cnt;
                bmax = cnt;
                bucket[cnt] = 1;
                for i in (cell1 + 1)..=cell2 {
                    cnt = (workset.clone() & g.0[lab[i]].clone()).count_ones();
                    while bmin > cnt {
                        bmin -= 1;
                        bucket[bmin] = 0;
                    }
                    while bmax < cnt {
                        bmax += 1;
                        bucket[bmax] = 0;
                    }
                    bucket[cnt] += 1;
                    count[i] = cnt;
                }
                if bmin == bmax
                {
                    longcode = mash(longcode,bmin+cell1);
                    continue;
                }
                c1 = cell1;
                maxcell = -1;
                for i in (bmin + 1)..=bmax {
                    if bucket[i] != 0 {
                        c2 = c1 + bucket[i];
                        bucket[i] = c1;
                        longcode = mash(longcode,i+c1);
                        if (c2 - c1) as isize > maxcell {
                            maxcell = (c2 - c1) as isize;
                            maxpos = Some(c1);
                        }
                        if c1 != cell1 {
                            active.add_one(c1);
                            if c2 - c1 == 1 {
                                hint = c1;
                            }
                            *numcells += 1;
                        }
                        if c2 <= cell2 {
                            ptn[c2 - 1] = level;
                            c1 = c2;
                        }
                    }
                }
                for i in (cell1+1)..=cell2 {
                    lab[i] = workperm[i];
                }
                if !active[cell1] {
                    active.add_one(cell1);
                    if let Some(maxpos) = maxpos {
                        active.remove_one(maxpos);
                    }
                }

            }}
        }
    }
    longcode = mash(longcode, *numcells);
    *code = cleanup(longcode);
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::nauty::{
        NAUTY_INFINITY, NAUTY_INFINITY_I,
        test::{create_diamond, create_n_circle, create_zero},
    };
    use bitvec::{bitvec, order::Msb0};

    #[test]
    fn test_is_autom() {
        assert!(is_autom(&create_n_circle(4), &[1, 2, 3, 0]));
        assert!(!is_autom(&create_n_circle(4), &[1, 0, 2, 3]));
    }

    #[test]
    fn test_refine_diamond_unpartitioned() {
        let mut g = create_diamond();
        let mut lab = [0, 1, 2, 3];
        let mut ptn = [1, 1, 1, 1, 0];
        let mut count = vec![4];
        let mut numcells = 1;
        let mut active = bitvec![usize, Msb0; 1; 4];
        let mut code: usize = 0;
        refine(
            &mut g,
            &mut lab,
            &mut ptn,
            0,
            &mut numcells,
            &mut count,
            &mut active,
            &mut code,
        );
    }

    #[test]
    fn test_3_0() {
        let mut g = create_zero(3);
        let mut lab = [0, 2, 1];
        let mut ptn = [2, NAUTY_INFINITY_I, 0];
        let mut count = vec![1, 0];
        let mut numcells = 2;
        let mut active = bitvec![usize, Msb0; 1, 0, 0];
        let mut code: usize = 21845;
        refine(
            &mut g,
            &mut lab,
            &mut ptn,
            2,
            &mut numcells,
            &mut count,
            &mut active,
            &mut code,
        );
        assert_eq!(lab, [0, 1, 2]);
        assert_eq!(ptn, [2, NAUTY_INFINITY_I, 0]);
        assert_eq!(count, [1, 0]);
        assert_eq!(numcells, 2);
        assert_eq!(active, bitvec![usize, Msb0; 0; 3]);
        assert_eq!(code, 4);
    }

    #[test]
    fn test_3_1() {
        let mut g = create_zero(3);
        let mut lab = [1, 0, 2];
        let mut ptn = [2, NAUTY_INFINITY_I, 0];
        let mut numcells = 2;
        let mut count = vec![0, 2];
        let mut active: bitvec::prelude::BitVec<usize, Msb0> = bitvec![usize, Msb0; 1, 0, 0];
        let mut code: usize = 1431812424;
        refine(
            &mut g,
            &mut lab,
            &mut ptn,
            2,
            &mut numcells,
            &mut count,
            &mut active,
            &mut code,
        );
        assert_eq!(lab, [1, 2, 0]);
        assert_eq!(ptn, [2, NAUTY_INFINITY_I, 0]);
        assert_eq!(count, [0, 2]);
        assert_eq!(numcells, 2);
        assert_eq!(active, bitvec![usize, Msb0; 0; 3]);
        assert_eq!(code, 4);
    }

    #[test]
    fn test_3_2() {
        let mut g = Graph::from_32(&[1610612736, 2684354560, 3221225472]);
        let mut lab = [0, 2, 1];
        let mut ptn = [2, NAUTY_INFINITY_I, 0];
        let mut numcells = 2;
        let mut count = vec![1, 0];
        let mut active: bitvec::prelude::BitVec<usize, Msb0> = bitvec![usize, Msb0; 1, 0, 0];
        let mut code: usize = 21845;
        refine(
            &mut g,
            &mut lab,
            &mut ptn,
            2,
            &mut numcells,
            &mut count,
            &mut active,
            &mut code,
        );
        assert_eq!(lab, [0, 2, 1]);
        assert_eq!(ptn, [2, NAUTY_INFINITY_I, 0]);
        assert_eq!(count, [1, 0]);
        assert_eq!(numcells, 2);
        assert_eq!(active, bitvec![usize, Msb0; 0; 3]);
        assert_eq!(code, 4);
    }

    #[test]
    fn test_3_3() {
        let mut g = Graph::from_32(&[1610612736, 2684354560, 3221225472]);
        let mut lab = [1, 0, 2];
        let mut ptn = [2, NAUTY_INFINITY_I, 0];
        let mut numcells = 2;
        let mut count = vec![0, 2];
        let mut active: bitvec::prelude::BitVec<usize, Msb0> = bitvec![usize, Msb0; 1, 0, 0];
        let mut code: usize = 1431812424;
        refine(
            &mut g,
            &mut lab,
            &mut ptn,
            2,
            &mut numcells,
            &mut count,
            &mut active,
            &mut code,
        );
        assert_eq!(lab, [1, 0, 2]);
        assert_eq!(ptn, [2, NAUTY_INFINITY_I, 0]);
        assert_eq!(count, [0, 2]);
        assert_eq!(numcells, 2);
        assert_eq!(active, bitvec![usize, Msb0; 0; 3]);
        assert_eq!(code, 4);
    }

    #[test]
    fn test_4_0() {
        let mut g = Graph::no_edge(4);
        let mut lab = [0, 3, 2, 1];
        let mut ptn = [2, NAUTY_INFINITY_I, NAUTY_INFINITY_I, 0];
        let mut numcells = 2;
        let mut count = vec![1, 0];
        let mut active: bitvec::prelude::BitVec<usize, Msb0> = bitvec![usize, Msb0; 1, 0, 0, 0];
        let mut code: usize = 21845;
        refine(
            &mut g,
            &mut lab,
            &mut ptn,
            2,
            &mut numcells,
            &mut count,
            &mut active,
            &mut code,
        );
        assert_eq!(lab, [0, 2, 1, 3]);
        assert_eq!(ptn, [2, NAUTY_INFINITY_I, NAUTY_INFINITY_I, 0]);
        assert_eq!(numcells, 2);
        assert_eq!(count, [1, 0]);
        assert_eq!(active, bitvec![usize, Msb0; 0; 4]);
        assert_eq!(code, 4);
    }
}
