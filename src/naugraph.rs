use crate::{
    nautil::SetWordNautilTrait,
    nauty::{
        Graph, Set, SetTrait,
        partition_nest::partition::{Partition, cell_mut::SplitResult},
    },
};
use bitvec::{bitvec, order::Msb0};
use cfor::cfor;
use itertools::Itertools;

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

/// Algorithm 1, p. 8 https://arxiv.org/pdf/1301.1493
///
/// while α is not empty and π is not discrete do
///     Remove some element W from α.
///     for each cell X of π do
///         Let X1 , . . . , Xk be the fragments of X distinguished according
///         to the number of edges from each vertex to W .
///         Replace X by X1 , . . . , Xk in π.
///         if X ∈ α then
///             Replace X by X1 , . . . , Xk in α.
///         else
///             Add all but one of the largest of X1 , . . . , Xk to α.
///         end
///     end
/// end
///
///    g: &mut Graph,   
///    lab: &mut [usize],        labels
///    ptn: &mut [isize],        partition
///    level: isize,             recursion level
///    numcells: &mut usize,     number of cells
///    count: &mut Vec<usize>,   number of vertices in cells
///    active: &mut Set,         vertices not fixed yet
///    code: &mut usize,         
fn refine_nest(
    g: &mut Graph,
    partition: &mut Partition,
    count: &mut [usize],
    active: &mut Set,
    code: &mut usize,
) {
    let mut i: usize;
    let mut split1: usize;
    let mut split2: usize;
    let mut cnt: usize;
    let mut bmin: usize;
    let mut bmax: usize;
    let mut longcode: usize;
    let mut maxcell: isize;
    let mut maxpos: Option<usize> = None;
    let mut hint: usize;
    let mut workperm: Vec<usize> = vec![0; g.n()];
    let mut bucket: Vec<usize> = vec![0; g.n() + 2];

    longcode = partition.numcells();
    hint = 0;
    loop {
        if partition.numcells() == g.n() {
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
        split2 = partition.cell_end(split1);
        let splitters = partition.get_splitters(split1);
        longcode = mash(longcode, split1 + split2);
        /* trivial splitting cell */
        if split1 == split2 {
            let gptr = &g.0[partition[split1]];
            let mut cells = partition.cells_mut();
            while let Some(mut cell) = cells.next() {
                if cell.len() == 1 {
                    continue;
                }
                let i = cell.split_trivial(gptr);
                if i > 0 && i < cell.len() {
                    longcode = mash(longcode, cell.partition_index(i - 1));
                    if active[cell.partition_index(i)] || 2 * i >= cell.len() {
                        active.add_one(cell.partition_index(i));
                        if i == cell.len() - 1 {
                            hint = cell.partition_index(i);
                        }
                    } else {
                        active.add_one(cell.partition_index(0));
                        if i == 1 {
                            hint = cell.partition_index(0);
                        }
                    }
                }
            }
        /* nontrivial splitting cell */
        } else {
            longcode = mash(longcode, split2 - split1 + 1);
            let mut cells = partition.cells_mut();
            while let Some(mut cell) = cells.next() {
                if cell.len() == 1 {
                    continue;
                }
                match cell.split_from_splitters(&splitters, g) {
                    SplitResult::Const(bmin) => {
                        longcode = mash(longcode, bmin + cell.partition_index(0));
                        continue;
                    }
                    SplitResult::Split {
                        biggest_cell_pos,
                        cells_indices,
                        cells_values,
                    } => {
                        for ((c1, c2), value) in
                            cells_indices.into_iter().tuple_windows().zip(cells_values)
                        {
                            longcode = mash(longcode, value + cell.partition_index(c1));
                            if c1 != 0 {
                                active.add_one(cell.partition_index(c1));
                                if c2 - c1 == 1 {
                                    hint = c1;
                                }
                            }
                        }
                        if !active[cell.partition_index(0)] {
                            active.add_one(cell.partition_index(0));
                            active.remove_one(cell.partition_index(biggest_cell_pos));
                        }
                    }
                }
            }
        }
    }
    longcode = mash(longcode, partition.numcells());
    *code = cleanup(longcode);
}

///
///
///    g: &mut Graph,   
///    lab: &mut [usize],        labels
///    ptn: &mut [isize],        partition
///    level: isize,             recursion level
///    numcells: &mut usize,     number of cells
///    count: &mut Vec<usize>,   number of vertices in cells
///    active: &mut Set,         vertices not fixed yet
///    code: &mut usize,         
#[allow(clippy::too_many_arguments)]
fn refine(
    g: &mut Graph,
    lab: &mut [usize],
    ptn: &mut [usize],
    level: usize,
    numcells: &mut usize,
    count: &mut Vec<usize>,
    active: &mut Set,
    code: &mut usize,
) {
    let mut i: usize;
    let mut c1: usize;
    let mut c2: isize;
    let mut labc1: usize;
    let mut split1: usize;
    let mut split2: usize;
    let mut cell1: usize;
    let mut cell2: usize;
    let mut cnt: usize;
    let mut bmin: usize;
    let mut bmax: usize;
    let mut longcode: usize;
    let mut maxcell: isize;
    let mut maxpos: Option<usize> = None;
    let mut hint: usize;
    let mut workperm: Vec<usize> = vec![0; g.n()];
    let mut bucket: Vec<usize> = vec![0; g.n() + 2];

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
            let gptr = &g.0[lab[split1]];
            cell2 = 0;
            cfor! {cell1 = 0; cell1 < g.n(); cell1 = cell2 + 1; {
                cfor! {cell2 = cell1; ptn[cell2] > level; cell2 += 1; {}}
                if cell1 == cell2 {
                    continue;
                }
                c1 = cell1;
                c2 = cell2 as isize;
                while c1 as isize <= c2 {
                    labc1 = lab[c1];
                    if gptr[labc1] {
                        c1 +=1;
                    } else {
                        lab[c1] = lab[c2 as usize];
                        lab[c2 as usize] = labc1;
                        c2 -= 1;
                    }
                }
                if c2 >= cell1 as isize && c1 <= cell2{
                    ptn[c2 as usize] = level;
                    longcode = mash(longcode, c2 as usize);
                    *numcells += 1;
                    if active[cell1] || (c2 as usize - cell1 >= cell2-c1) {
                        active.add_one(c1);
                        if c1 == cell2 {
                            hint = c1;
                        }
                    } else {
                        active.add_one(cell1);
                        if c2 as usize == cell1 {
                            hint = cell1;
                        }
                    }
                }
            }}
        /* nontrivial splitting cell */
        } else {
            let mut workset: Set = bitvec![usize, Msb0; 0; g.n()];
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
                let mut c2: usize;
                for i in bmin..=bmax {
                    if bucket[i] != 0 {
                        c2 = c1 + bucket[i];
                        bucket[i] = c1;
                        longcode = mash(longcode, i + c1);
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
                for i in cell1..=cell2 {
                    workperm[bucket[count[i]]] = lab[i];
                    bucket[count[i]] += 1;
                }
                lab[(cell1 + 1)..(cell2 + 1)].copy_from_slice(&workperm[(cell1 + 1)..(cell2 + 1)]);
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
        NAUTY_INFINITY,
        partition_nest::PartitionNest,
        test::{
            create_complete, create_diamond, create_from_offsets, create_pentagram, graphs::cycle,
        },
        u32_to_bitvec,
    };
    use bitvec::{bitvec, order::Msb0, vec::BitVec};
    use test_case::test_case;

    #[test]
    fn test_is_autom() {
        assert!(is_autom(&cycle(4), &[1, 2, 3, 0]));
        assert!(!is_autom(&cycle(4), &[1, 0, 2, 3]));
    }

    //  test_no_nest
    // #[test_case(create_diamond(), &[0, 1, 2, 3], &[NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0], 0, 1, &[4, 4, 4, 4], bitvec![usize, Msb0; 1; 4], 0, &[0, 2, 3, 1], &[0, 0, NAUTY_INFINITY, 0], 2, &[3, 2, 0, 0], bitvec![usize, Msb0; 0; 4], 27493; "diamond_unpartitioned")]
    #[test_case(create_from_offsets(7, &[3, 4]),  &[0,3,4,2,1,6,5], &[2, NAUTY_INFINITY, 2, NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, 3, &[5, 0, 0, 0, 1, 1, 0], bitvec![usize, Msb0; 0,1,0,0,0,0,0], 21845, &[0,3,4,2,5,1,6], &[2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 0], 4, &[5,1,1,1,1,0,0], bitvec![usize, Msb0; 0; 7], 27483; "FCp`_")]
    #[test_case(create_pentagram(),  &[1, 0, 3, 2, 4], &[2, NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, 2, &[0, 4, 3, 2, 1], u32_to_bitvec(2147483648, 5), 1431812424, &[1, 4, 3, 2, 0], &[2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 0], 3, &[0,1, 1, 1, 1], bitvec![usize, Msb0; 0; 5], 27427; "5_2")]
    #[test_case(create_pentagram(),  &[0, 4, 3, 2, 1], &[2, NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, 2, &[3, 2, 1, 0,4], u32_to_bitvec(2147483648, 5), 21845, &[0, 2, 3, 1, 4], &[2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 0], 3, &[3,1,1,1,1], bitvec![usize, Msb0; 0; 5], 27427; "DUW 1")]
    #[test_case(cycle(3), &[0, 2, 1], &[2, NAUTY_INFINITY, 0], 2, 2, &[1, 0], bitvec![usize, Msb0; 1, 0, 0], 1431812424, &[0, 2, 1], &[2, NAUTY_INFINITY, 0], 2, &[1, 0], bitvec![usize, Msb0; 0; 3], 4; "3_2")]
    #[test_case(cycle(3), &[1, 0, 2], &[2, NAUTY_INFINITY, 0], 2, 2, &[0, 2], bitvec![usize, Msb0; 1, 0, 0], 1431812424, &[1, 0, 2], &[2, NAUTY_INFINITY, 0], 2, &[0, 2], bitvec![usize, Msb0; 0; 3], 4; "3_3")]
    #[test_case(Graph::no_edge(3), &[1, 0, 2], &[2, NAUTY_INFINITY, 0], 2, 2, &[0, 2], bitvec![usize, Msb0; 1, 0, 0], 1431812424, &[1, 2, 0], &[2, NAUTY_INFINITY, 0], 2, &[0, 2], bitvec![usize, Msb0; 0; 3], 4; "3_1")]
    #[test_case(Graph::no_edge(3), &[0, 2, 1], &[2, NAUTY_INFINITY, 0], 2, 2, &[1, 0], bitvec![usize, Msb0; 1, 0, 0], 21845, &[0, 1, 2], &[2, NAUTY_INFINITY, 0], 2, &[1, 0], bitvec![usize, Msb0; 0; 3], 4)]
    #[test_case(Graph::no_edge(4), &[0, 3, 2, 1], &[2, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, 2, &[1, 0], bitvec![usize, Msb0; 1, 0, 0, 0], 21845, &[0, 2, 1, 3], &[2, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, &[1, 0], bitvec![usize, Msb0; 0; 4], 4)]
    #[test_case(create_pentagram(), &[0, 2, 3, 1, 4], &[3, 3, 2, NAUTY_INFINITY, 0], 3, 4, &[3, 1, 1], u32_to_bitvec(1073741824, 5), 21845, &[0, 2, 3, 4, 1], &[3, 3, 2, 3, 0], 5, &[3, 1, 1], u32_to_bitvec(134217728, 5), 27417; "5_0")]
    #[test_case(create_complete(4), &[1, 0, 2, 3], &[2, 3, NAUTY_INFINITY, 0], 3, 3, &[0, 2, 1], u32_to_bitvec(1073741824, 4), 1431812424, &[1, 0, 2, 3], &[2, 3, NAUTY_INFINITY, 0], 3, &[0, 2, 1], bitvec![usize, Msb0; 0; 4], 64)]
    #[test_case(Graph::no_edge(4), &[1, 0, 2, 3], &[2, NAUTY_INFINITY, NAUTY_INFINITY, 0], 3, 2, &[0, 2], u32_to_bitvec(2147483648, 4), 1431812424, &[1, 2, 3, 0], &[2, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, &[0, 2], bitvec![usize, Msb0; 0; 4], 4; "4_3")]
    #[test_case(Graph::no_edge(4), &[0, 2, 1, 3], &[2, 3, NAUTY_INFINITY, 0], 3, 3, &[0, 1, 3], u32_to_bitvec(1073741824, 4), 1431812424, &[0, 2, 3, 1], &[2, 3, NAUTY_INFINITY, 0], 3, &[0, 1, 3], bitvec![usize, Msb0; 0; 4], 64)]
    #[test_case(Graph::no_edge(4), &[0, 1, 2, 3], &[2, 3, NAUTY_INFINITY, 0], 3, 3, &[1, 0, 2], u32_to_bitvec(1073741824, 4), 21845, &[0, 1, 3, 2], &[2, 3, NAUTY_INFINITY, 0], 3, &[1, 0, 2], bitvec![usize, Msb0; 0; 4], 64)]
    #[test_case(Graph::no_edge(4), &[1, 0, 2, 3], &[2, 3, NAUTY_INFINITY, 0], 3, 3, &[0, 2, 1], u32_to_bitvec(1073741824, 4), 1431812424, &[1, 0, 3, 2], &[2, 3, NAUTY_INFINITY, 0], 3, &[0, 2, 1], bitvec![usize, Msb0; 0; 4], 64)]
    #[allow(clippy::too_many_arguments)]
    fn test_no_nest(
        mut g: Graph,
        lab: &[usize],
        ptn: &[usize],
        level: usize,
        expected_numcells_before: usize,
        count: &[usize],
        mut active: BitVec<usize, Msb0>,
        mut code: usize,
        expected_lab: &[usize],
        expected_ptn: &[usize],
        expected_numcells_after: usize,
        expected_count: &[usize],
        expected_active: BitVec<usize, Msb0>,
        expected_code: usize,
    ) {
        println!("{}", g.to_matrix());
        println!("{}", g.to_graph6());
        let mut lab = lab.to_vec();
        let mut ptn = ptn.to_vec();
        assert_eq!(
            Partition::new(PartitionNest::new(lab.clone(), ptn.clone()), level).numcells(),
            expected_numcells_before
        );
        let mut count = count.to_vec();
        let mut numcells = expected_numcells_before;
        refine(
            &mut g,
            &mut lab,
            &mut ptn,
            level,
            &mut numcells,
            &mut count,
            &mut active,
            &mut code,
        );
        assert_eq!(lab, expected_lab);
        assert_eq!(ptn, expected_ptn);
        assert_eq!(numcells, expected_numcells_after);
        assert_eq!(count, expected_count);
        assert_eq!(active, expected_active);
        assert_eq!(code, expected_code);
    }

    //  test_nest
    //#[test_case(create_diamond(), &[0, 1, 2, 3], &[NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0], 0, 1, &[4, 4, 4, 4], bitvec![usize, Msb0; 1; 4], 0, &[0, 2, 3, 1], &[0, 0, NAUTY_INFINITY, 0], 2, &[3, 2, 0, 0], bitvec![usize, Msb0; 0; 4], 27493; "diamond_unpartitioned")]
    #[test_case(create_from_offsets(7, &[3, 4]), &[0,3,4,2,1,6,5], &[2, NAUTY_INFINITY, 2, NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, 3, &[5, 0, 0, 0, 1, 1, 0], bitvec![usize, Msb0; 0,1,0,0,0,0,0], 21845, &[0,3,4,2,5,1,6], &[2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 0], 4,  bitvec![usize, Msb0; 0; 7], 27483; "FCp`_")]
    #[test_case(create_pentagram(),  &[1, 0, 3, 2, 4], &[2, NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, 2, &[0, 4, 3, 2, 1], u32_to_bitvec(2147483648, 5), 1431812424, &[1, 4, 3, 2, 0], &[2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 0], 3,  bitvec![usize, Msb0; 0; 5], 27427; "DUW 2")]
    #[test_case(create_pentagram(),  &[0, 4, 3, 2, 1], &[2, NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, 2, &[3, 2, 1, 0,4], u32_to_bitvec(2147483648, 5), 21845, &[0, 2, 3, 1, 4], &[2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 0], 3, bitvec![usize, Msb0; 0; 5], 27427; "DUW 1")]
    #[test_case(cycle(3), &[0, 2, 1], &[2, NAUTY_INFINITY, 0], 2, 2, &[1, 0], bitvec![usize, Msb0; 1, 0, 0], 1431812424, &[0, 2, 1], &[2, NAUTY_INFINITY, 0], 2,  bitvec![usize, Msb0; 0; 3], 4; "3_2")]
    #[test_case(cycle(3), &[1, 0, 2], &[2, NAUTY_INFINITY, 0], 2, 2, &[0, 2], bitvec![usize, Msb0; 1, 0, 0], 1431812424, &[1, 0, 2], &[2, NAUTY_INFINITY, 0], 2,  bitvec![usize, Msb0; 0; 3], 4; "3_3")]
    #[test_case(Graph::no_edge(3), &[1, 0, 2], &[2, NAUTY_INFINITY, 0], 2, 2, &[0, 2], bitvec![usize, Msb0; 1, 0, 0], 1431812424, &[1, 2, 0], &[2, NAUTY_INFINITY, 0], 2,  bitvec![usize, Msb0; 0; 3], 4; "3_1")]
    #[test_case(Graph::no_edge(3), &[0, 2, 1], &[2, NAUTY_INFINITY, 0], 2, 2, &[1, 0], bitvec![usize, Msb0; 1, 0, 0], 21845, &[0, 1, 2], &[2, NAUTY_INFINITY, 0], 2,  bitvec![usize, Msb0; 0; 3], 4)]
    #[test_case(Graph::no_edge(4), &[0, 3, 2, 1], &[2, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2, 2, &[1, 0], bitvec![usize, Msb0; 1, 0, 0, 0], 21845, &[0, 2, 1, 3], &[2, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2,  bitvec![usize, Msb0; 0; 4], 4)]
    #[test_case(create_pentagram(), &[0, 2, 3, 1, 4], &[3, 3, 2, NAUTY_INFINITY, 0], 3, 4, &[3, 1, 1], u32_to_bitvec(1073741824, 5), 21845, &[0, 2, 3, 4, 1], &[3, 3, 2, 3, 0], 5,  u32_to_bitvec(134217728, 5), 27417; "5_0")]
    #[test_case(Graph::from_u32(&[1879048192, 2952790016, 3489660928, 3758096384]), &[1, 0, 2, 3], &[2, 3, NAUTY_INFINITY, 0], 3, 3, &[0, 2, 1], u32_to_bitvec(1073741824, 4), 1431812424, &[1, 0, 2, 3], &[2, 3, NAUTY_INFINITY, 0], 3,  bitvec![usize, Msb0; 0; 4], 64)]
    #[test_case(Graph::no_edge(4), &[1, 0, 2, 3], &[2, 3, NAUTY_INFINITY, 0], 3, 3, &[0, 2, 1], u32_to_bitvec(1073741824, 4), 1431812424, &[1, 0, 3, 2], &[2, 3, NAUTY_INFINITY, 0], 3,  bitvec![usize, Msb0; 0; 4], 64)]
    #[test_case(Graph::no_edge(4), &[1, 0, 2, 3], &[2, NAUTY_INFINITY, NAUTY_INFINITY, 0], 3, 2, &[0, 2], u32_to_bitvec(2147483648, 4), 1431812424, &[1, 2, 3, 0], &[2, NAUTY_INFINITY, NAUTY_INFINITY, 0], 2,  bitvec![usize, Msb0; 0; 4], 4; "4_3")]
    #[test_case(Graph::no_edge(4), &[0, 2, 1, 3], &[2, 3, NAUTY_INFINITY, 0], 3, 3, &[0, 1, 3], u32_to_bitvec(1073741824, 4), 1431812424, &[0, 2, 3, 1], &[2, 3, NAUTY_INFINITY, 0], 3,  bitvec![usize, Msb0; 0; 4], 64)]
    #[test_case(Graph::no_edge(4), &[0, 1, 2, 3], &[2, 3, NAUTY_INFINITY, 0], 3, 3, &[1, 0, 2], u32_to_bitvec(1073741824, 4), 21845, &[0, 1, 3, 2], &[2, 3, NAUTY_INFINITY, 0], 3,  bitvec![usize, Msb0; 0; 4], 64)]
    #[allow(clippy::too_many_arguments)]
    fn test_nest(
        mut g: Graph,
        lab: &[usize],
        ptn: &[usize],
        level: usize,
        expected_numcells_before: usize,
        count: &[usize],
        mut active: BitVec<usize, Msb0>,
        mut code: usize,
        expected_lab: &[usize],
        expected_ptn: &[usize],
        expected_numcells_after: usize,
        expected_active: BitVec<usize, Msb0>,
        expected_code: usize,
    ) {
        println!("{}", g.to_matrix());
        println!("{}", g.to_graph6());
        let nest = PartitionNest::new(lab.to_vec(), ptn.to_vec());
        let mut partition = Partition::new(nest, level);
        assert_eq!(partition.numcells(), expected_numcells_before);
        let mut count = count.to_vec();
        refine_nest(&mut g, &mut partition, &mut count, &mut active, &mut code);
        let expected_nest = PartitionNest::new(expected_lab.to_vec(), expected_ptn.to_vec());
        let expected_partition = Partition::new(expected_nest, level);
        assert_eq!(partition, expected_partition);
        assert_eq!(partition.numcells(), expected_numcells_after);
        //assert_eq!(count, expected_count);
        assert_eq!(active, expected_active);
        assert_eq!(code, expected_code);
    }
}
