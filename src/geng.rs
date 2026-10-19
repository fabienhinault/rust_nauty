use crate::{
    cli::GengCli,
    graph::Graph,
    level_data::LevelData,
    nauty::{
        NautyCounter, Set, SetTrait, graph,
        partition_nest::{
            PartitionNest,
            partition::{Partition, cell_mut::SplitResult},
        },
    },
};
use bitvec::prelude::Msb0;
use bitvec::{prelude::BitVec, view::BitView};
use itertools::Itertools;

#[inline(always)]
fn xbit(i: usize) -> usize {
    1 << i
}

fn xnextbit(i: usize) -> usize {
    i.leading_zeros() as usize
}

impl Graph {
    /* Return TRUE iff there is a K4 including the last vertex */
    // geng.c 1469
    pub fn hask4(&self) -> bool {
        todo!();
    }

    /* Return TRUE if there is a claw (induced K(1,3)) involving the last vertex */
    // geng.c 1492
    pub fn hasclaw(&self) -> bool {
        todo!();
    }

    /* This function speeds up the generation of connected graphs
    with not many edges. */
    pub fn connpreprune(&self, args: &GengCli) -> bool {
        let maxn = args.maxn;
        let maxe = args.maxe();
        if self.n() == maxn || maxe - maxn >= 5 {
            return false;
        }

        let ne = self.edge_count();
        let mut nc = 0;
        let mut notvisited = Set::ones(self.n());
        while notvisited.any() {
            nc += 1;
            let mut queue = notvisited.masked_from(notvisited.last_one().expect("last_one"));
            while queue.any() {
                let i = queue.takebit();
                notvisited.set(i, false);
                queue |= self.0[i].clone() & notvisited.clone();
            }
        }
        ne - self.n() + nc > maxe - maxn + 1
    }
}

pub struct Extender {
    pub maxn: usize,
    pub mindeg: usize,
    pub connec: u8,
    pub nodes: Vec<NautyCounter>,
    pub data: Vec<LevelData>,
    pub canonise: bool,
    pub args: GengCli,
}

impl Extender {
    /* decide if n in theta(g+x) -  version for n+1 < maxn */
    // 1711
    fn accept1(
        &self,
        g: &mut graph::Graph,
        x: usize,
        gx: &mut graph::Graph,
        deg: &mut [usize],
        d: bool,
    ) -> bool {
        /* decide if n in theta(g+x) -  version for n+1 < maxn */
        let n = g.n();
        let nx = n + 1;
        gx.0[..n].clone_from_slice(&g.0[..n]);
        gx.0[n] = Set::zeros(nx);
        let degn = x.count_ones() as usize;
        deg[n] = degn;
        let mut xw = BitVec::from_bitslice(x.view_bits::<Msb0>().split_at(n).0);
        while xw.any() {
            let i = xw.leading_zeros() as usize;
            let xwi = xw[i];
            xw.set(i, !xwi);
            gx[i].set(n, true);
            gx[n].set(i, true);
            deg[i] += 1;
        }
        if self.args.k4free && gx.hask4() {
            return false;
        }
        if self.args.clawfree && gx.hasclaw() {
            return false;
        }
        if self.args.connec2 && n + 2 == self.maxn && !g.isconnected() {
            return false;
        }
        if ((self.args.connec2 && n + 2 < self.maxn) || (self.args.connec1 && n + 2 <= self.maxn))
            && gx.connpreprune(&self.args)
        {
            return false;
        }
        let nest = PartitionNest::from_deg(deg);
        let mut active = Set::zeros(nx);
        active.set(0, true);
        if let Some(&extra_zero_index) = nest.first_extra_zero() {
            active.set(extra_zero_index + 1, true);
        }

        false
    }

    fn accept2(
        &self,
        g: &mut graph::Graph,
        n: usize,
        x: usize,
        gx: &mut graph::Graph,
        deg: &[usize],
        d: bool,
    ) -> bool {
        false
    }
    /* extend from n to n+1 -- version for general graphs */
    // 2314
    pub fn genextend(
        &mut self,
        g: &mut graph::Graph,
        n: usize,
        deg: &mut [usize],
        ne: usize,
        rigid: bool,
        xlb: usize,
        xub: usize,
        ecount: &mut [i32],
    ) {
        let mut xlb = xlb;

        self.nodes[n] += 1;
        let nx = n + 1;
        let dmax = deg[n - 1];
        let dcrit = self.mindeg - self.maxn + n;
        let mut d: usize = 0;
        let mut dlow = 0;
        for i in 0..n {
            if deg[i] == dmax {
                d |= xbit(i);
            }
            if deg[i] == dcrit {
                dlow |= xbit(i);
            }
        }
        if xlb == dmax && d.count_ones() as usize + dmax > n {
            xlb += 1;
        }
        if nx == self.maxn && xlb < self.mindeg {
            xlb = self.mindeg;
        }
        if xlb > xub {
            return;
        }

        let imin = self.data[n].xstart[xlb];
        let imax = self.data[n].xstart[xub + 1];
        //let x_set_card = &mut data.x_set_card;
        let xorb = &self.data[n].xorb;
        let mut gx: graph::Graph = graph::Graph::default();
        if nx == self.maxn {
            for i in imin..imax {
                if !rigid && xorb[i] != i {
                    continue;
                }
                let x = self.data[n].xset(i);
                let xc = self.data[n].xcard(i);
                if xc == dmax && (x & d) != 0 {
                    continue;
                }
                if (dlow & !x) != 0 {
                    continue;
                }
                if self.accept2(
                    g,
                    n,
                    x,
                    &mut gx,
                    deg,
                    xc > dmax + 1 || (xc == dmax + 1 && (x & d) == 0),
                ) && (self.connec == 0
                    || (self.connec == 1 && gx.isconnected())
                    || (self.connec == 2 && gx.isbiconnected()))
                {
                    ecount[ne + xc] += 1;
                    println!(
                        "{}",
                        if self.canonise {
                            gx.canonise().to_graph6()
                        } else {
                            gx.to_graph6()
                        }
                    );
                }
            }
        } else {
            for i in imin..imax {
                if !rigid && xorb[i] != i {
                    continue;
                }
                let x = self.data[n].xset(i);
                let xc = self.data[n].xcard(i);
                if xc == dmax && (x & d) != 0 {
                    continue;
                }
                if (dlow & !x) != 0 {
                    continue;
                }
                if self.accept1(
                    g,
                    x,
                    &mut gx,
                    deg,
                    xc > dmax + 1 || (xc == dmax + 1 && (x & d) == 0),
                ) && (self.connec == 0
                    || (self.connec == 1 && gx.isconnected())
                    || (self.connec == 2 && gx.isbiconnected()))
                {
                    ecount[ne + xc] += 1;
                    println!(
                        "{}",
                        if self.canonise {
                            gx.canonise().to_graph6()
                        } else {
                            gx.to_graph6()
                        }
                    );
                }
            }
        }
    }
}

/* form orbits on powerset of VG
called by nauty;  operates on data[n] */
// 1190
fn userautomproc(
    count: usize,
    p: &[usize],
    orbits: &[usize],
    numorbits: usize,
    stabvertex: usize,
    n: usize,
    data: &mut LevelData,
) {
    // let xorb = &mut data.xorb;
    // let xset = |i| data.xset(i);
    // let xinv = &mut data.xinv;
    let LevelData {
        xorb,
        x_set_card,
        xinv,
        ..
    } = data;
    let xset = |i: usize| x_set_card[i].set;
    let lo = data.lo;
    let hi = data.hi;

    if count == 1 {
        for i in lo..hi {
            xorb[i] = i;
        }
    }
    let mut moved = 0;
    for i in 0..n {
        if p[i] != i {
            moved |= 1 << i;
        }
    }
    for i in lo..hi {
        let mut w = xset(i);
        if w & moved == 0 {
            continue;
        }
        let mut pxi = xset(i) & !moved;
        while w != 0 {
            let j1 = xnextbit(w);
            w ^= 1 << j1;
            pxi |= 1 << p[j1];
        }
        let pi = xinv[pxi].expect("xinv[pxi]");
        let mut j1 = xorb[i];
        while xorb[j1] != j1 {
            j1 = xorb[j1];
        }
        let mut j2 = xorb[pi];
        while xorb[j2] != j2 {
            j2 = xorb[j2];
        }

        if j1 < j2 {
            xorb[pi] = j1;
            xorb[i] = j1;
            xorb[j2] = j1;
        } else if j1 > j2 {
            xorb[pi] = j2;
            xorb[i] = j2;
            xorb[j1] = j2;
        }
    }
}

/*****************************************************************************
*                                                                            *
*  refinex(g,lab,ptn,level,numcells,count,active,goodret,code,m,n) is a      *
*  custom version of refine() which can exit quickly if required.            *
*                                                                            *
*  Only use at level==0.                                                     *
*  goodret : whether to do an early return for code 1                        *
*  code := -1 for n-1 not max, 0 for maybe, 1 for definite                   *
*                                                                            *
*****************************************************************************/
pub fn refinex(
    g: &graph::Graph,
    partition: &mut Partition,
    active: &mut Set,
    goodret: bool,
    code: &mut isize,
) {
    let mut split1: usize;
    let mut split2: usize;
    let mut longcode: usize;
    let mut hint: usize;

    let n = g.n();
    if n == 1 {
        *code = 1;
        return;
    }
    *code = 0;
    let mut lact = active.clone();
    while !partition.is_discrete() && !lact.is_zero() {
        let split1 = lact.takebit();
        split2 = partition.cell_end(split1);
        let splitters = partition.get_splitters(split1);
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
                    lact.set(i, true);
                }
            }
        /* nontrivial splitting cell */
        } else {
            let mut cells = partition.cells_mut();
            while let Some(mut cell) = cells.next() {
                if cell.len() == 1 {
                    continue;
                }
                match cell.split_from_splitters(&splitters, g) {
                    SplitResult::Const(bmin) => {
                        continue;
                    }
                    SplitResult::Split {
                        cells_indices,
                        cells_values,
                        ..
                    } => {
                        for ((c1, _), _) in
                            cells_indices.into_iter().tuple_windows().zip(cells_values)
                        {
                            if c1 != 0 {
                                lact.add_one(cell.partition_index(c1));
                            }
                        }
                    }
                }
            }
        }
        if partition.nest.get_ptn(n - 2) == 0 {
            if partition[n - 1] == n - 1 {
                *code = 1;
                if goodret {
                    return;
                }
            } else {
                *code = -1;
                return;
            }
        } else {
            let mut i = n - 1;
            loop {
                if partition[i] == n - 1 {
                    break;
                }
                i -= 1;
                if partition.nest.get_ptn(i) == 0 {
                    *code = -1;
                    return;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
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

    #[test_case(Graph::from_u32(&[33554432,0,0,0,0,0,2147483648]), &[1,2,3,4,5,0,6], &[1,1,1,1,0,1,0], 0, 2, u32_to_bitvec(2214592512, 7), true, 0, &[1,2,3,4,5,0,6], &[1,1,1,1,0,1,0], 2, 0; "F??C?") ]
    #[test_case(Graph::no_edge(3), &[2,1,0], &[1,1,0], 0, 1, u32_to_bitvec(2147483648, 3), false, 1, &[2,1,0], &[1,1,0], 1, 0)]
    #[test_case(Graph::no_edge(2), &[1,0], &[1, 0], 0, 1, u32_to_bitvec(2147483648, 2), false, 32767, &[1, 0], &[1,0], 1, 0)]
    #[allow(clippy::too_many_arguments)]
    fn test_nest(
        g: graph::Graph,
        lab: &[usize],
        ptn: &[usize],
        level: usize,
        expected_numcells_before: usize,
        //count: &[usize],
        mut active: BitVec<usize, Msb0>,
        goodret: bool,
        mut code: isize,
        expected_lab: &[usize],
        expected_ptn: &[usize],
        expected_numcells_after: usize,
        expected_code: isize,
    ) {
        println!("{}", g.to_matrix());
        println!("{}", g.to_graph6());
        let nest = PartitionNest::new(lab.to_vec(), ptn.to_vec());
        let mut partition = Partition::new(nest, level);
        assert_eq!(partition.numcells(), expected_numcells_before);
        refinex(&g, &mut partition, &mut active, goodret, &mut code);
        let expected_nest = PartitionNest::new(expected_lab.to_vec(), expected_ptn.to_vec());
        let expected_partition = Partition::new(expected_nest, level);
        assert_eq!(partition, expected_partition);
        assert_eq!(partition.numcells(), expected_numcells_after);
        assert_eq!(code, expected_code);
    }
}
