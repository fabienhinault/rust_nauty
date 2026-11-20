use std::assert_eq;

use crate::{
    cli::GengCli,
    env::Env,
    graph::Graph,
    level_data::{LevelData, xnextbit},
    nauty::{
        NautyCounter, OptionBlk, Set, SetTrait, StatBlk, graph, nauty,
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

impl Graph {
    /* Return TRUE iff there is a K4 including the last vertex */
    // geng.c 1469
    pub fn hask4(&self) -> bool {
        todo!();
        false
    }

    /* Return TRUE if there is a claw (induced K(1,3)) involving the last vertex */
    // geng.c 1492
    pub fn hasclaw(&self) -> bool {
        todo!();
        false
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
    pub canonise: bool,
    pub args: GengCli,
    pub splitlevel: isize,
    pub env: Env,
}

impl Extender {
    /* decide if n in theta(g+x) -  version for n+1 < maxn */
    // 1711
    fn accept1(&self, g: &Graph, x: usize, data: &mut LevelData) -> Option<(Graph, bool)> {
        let n = g.n();
        let nx = n + 1;
        let mut gx = g.clone();
        gx.push(Set::from_sym_usize(x, nx));
        if self.args.k4free && gx.hask4() {
            return None;
        }
        if self.args.clawfree && gx.hasclaw() {
            return None;
        }
        if self.args.connec2 && n + 2 == self.maxn && !g.isconnected() {
            return None;
        }
        if ((self.args.connec2 && n + 2 < self.maxn) || (self.args.connec1 && n + 2 <= self.maxn))
            && gx.connpreprune(&self.args)
        {
            return None;
        }
        let nest = gx.degs_partition_nest();
        let mut active = Set::zeros(nx);
        active.set(0, true);
        if let Some(extra_zero_index) = nest.first_extra_zero() {
            active.set(extra_zero_index + 1, true);
        }
        let mut code = 0;
        let mut partition = Partition { nest, level: 0 };
        refinex(&gx, &mut partition, &mut active, false, &mut code);
        if code < 0 {
            return None;
        }
        if partition.is_discrete() {
            return Some((gx, true));
        }
        let mut options: OptionBlk = OptionBlk {
            getcanon: 1,
            defaultptn: false,
            userautomproc: Some(data),
            tc_level: 100,
            ..Default::default()
        };
        active = BitVec::zeros(nx);
        let mut orbits: Vec<usize> = vec![];
        let mut h: Graph = Graph::no_edge(nx);
        let mut stats: StatBlk = StatBlk::default();
        let _ = nauty(
            &mut gx,
            &mut partition.nest,
            Some(&active),
            &mut orbits,
            &mut options,
            &mut stats,
            &mut h,
        );
        options.userautomproc = None;
        if orbits[partition[n]] == orbits[n] {
            Some((gx, stats.numorbits == nx))
        } else {
            None
        }
    }

    /* decide if n in theta(g+x)  --  version for n+1 == maxn */
    // 1935
    #[allow(clippy::too_many_arguments)]
    fn accept2(
        &self,
        g: &Graph,
        x: usize,
        nuniq: bool,
        data: &mut LevelData,
    ) -> Option<(Graph, Option<Graph>)> {
        let n = g.n();
        let nx = n + 1;
        let mut gx = g.clone();
        gx.push(Set::from_sym_usize(x, nx));
        if self.args.k4free && gx.hask4() {
            return None;
        }
        if self.args.clawfree && gx.hasclaw() {
            return None;
        }
        if self.args.connec2 && n + 2 == self.maxn && !g.isconnected() {
            return None;
        }
        if ((self.args.connec2 && n + 2 < self.maxn) || (self.args.connec1 && n + 2 <= self.maxn))
            && gx.connpreprune(&self.args)
        {
            return None;
        }
        if nuniq {
            let gcan = self.args.canonise.then(|| gx.canonise());
            return Some((gx, gcan));
        }
        let mut nest = gx.degs_partition_nest();
        let mut active = Set::zeros(nx);
        active.set(0, true);
        match nest.first_extra_zero() {
            None => {
                let invar = gx.distinvar();
                match invar {
                    None => return None,
                    Some(invar) => {
                        let j0 = nest.apply_distance_invariant(0, &invar);
                        if j0 > 0 {
                            if j0 == n {
                                let gcan = self.args.canonise.then(|| gx.canonise());
                                return Some((gx, gcan));
                            }
                            active.set(j0, true);
                        }
                    }
                }
            }
            Some(extra_zero_index) => {
                active.set(extra_zero_index + 1, true);
                let mut vmax: Set = Set::zeros(nx);
                for i in (extra_zero_index + 1)..nx {
                    vmax.set(nest[i], true);
                }
                let mut gvs = vec![0; extra_zero_index + 1];
                gvs.append(
                    &mut ((extra_zero_index + 1)..nx)
                        .map(|i| (gx[i].clone() & vmax.clone()).count_ones())
                        .collect(),
                );

                let j0 = nest.apply_distance_invariant(extra_zero_index + 1, &gvs);
                if j0 > 0 {
                    if j0 == n {
                        let gcan = self.args.canonise.then(|| gx.canonise());
                        return Some((gx, gcan));
                    }
                    active.set(j0, true);
                }
            }
        }
        let mut code = 0;
        let mut partition = Partition { nest, level: 0 };
        refinex(&mut gx, &mut partition, &mut active, false, &mut code);
        if code < 0 {
            return None;
        }
        let cheapacc: bool = partition.cheapacc(code);
        if cheapacc {
            let gcan = self.args.canonise.then(|| gx.canonise());
            return Some((gx, gcan));
        }

        let mut options: OptionBlk = OptionBlk::default();
        options.getcanon = 1;
        options.defaultptn = false;
        options.userautomproc = Some(data);
        active = BitVec::zeros(nx);
        let mut orbits: Vec<usize> = vec![];
        let mut h: Graph = Graph::no_edge(nx);
        let mut stats: StatBlk = StatBlk::default();
        let _ = nauty(
            &mut gx,
            &mut partition.nest,
            Some(&active),
            &mut orbits,
            &mut options,
            &mut stats,
            &mut h,
        );
        options.userautomproc = None;
        if orbits[partition[n]] == orbits[n] {
            let gcan = self.args.canonise.then(|| gx.canonise());
            Some((gx, gcan))
        } else {
            None
        }
    }

    /* extend from n to n+1 -- version for general graphs */
    // 2314
    #[allow(clippy::too_many_arguments)]
    pub fn genextend(
        &mut self,
        g: &mut graph::Graph,
        ne: usize,
        rigid: bool,
        xlb: usize,
        xub: usize,
        ecount: &mut [usize],
        data: &mut Vec<LevelData>,
    ) {
        let (n, nx, dmax, d, dlow, imin, imax, xorb) = match self.init_extend(g, xlb, xub, data) {
            Some(value) => value,
            None => return,
        };
        if nx == self.maxn {
            for i in imin..imax {
                if !rigid && xorb[i] != i {
                    continue;
                }
                let x = data[n].xset(i);
                let xc = data[n].xcard(i);
                if xc == dmax && (x & d) != 0 {
                    continue;
                }
                if (dlow & !x) != 0 {
                    continue;
                }
                if let Some((gx, _)) = self.accept2(
                    g,
                    x,
                    xc > dmax + 1 || (xc == dmax + 1 && (x & d) == 0),
                    &mut data[n],
                ) && (self.connec == 0
                    || (self.connec == 1 && gx.isconnected())
                    || (self.connec == 2 && gx.isbiconnected()))
                {
                    ecount[ne + xc] += 1;
                    print!(
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
            let xs: Vec<usize> = (imin..imax).map(|i| data[n].xset(i)).collect_vec();
            let xcs: Vec<usize> = (imin..imax).map(|i| data[n].xcard(i)).collect_vec();
            for i in imin..imax {
                if !rigid && xorb[i] != i {
                    continue;
                }
                let x = xs[i - imin];
                let xc = xcs[i - imin];
                if xc == dmax && (x & d) != 0 {
                    continue;
                }
                if (dlow & !x) != 0 {
                    continue;
                }
                if nx as isize == self.splitlevel {
                    todo!();
                }
                if data[nx].edge_count != ne + xc || data[nx].dmax != xc {
                    data[nx].xbnds(&self.env, nx, ne + xc, xc);
                }
                let xlbx = data[nx].xlb;
                let xubx = data[nx].xub;
                if xlbx > xubx {
                    continue;
                }
                data[nx].lo = data[nx].xstart[xlbx];
                data[nx].hi = data[nx].xstart[xubx + 1];
                if let Some((mut gx, rigidx)) = self.accept1(g, x, &mut data[nx])
                    && (self.connec == 0
                        || (self.connec == 1 && gx.isconnected())
                        || (self.connec == 2 && gx.isbiconnected()))
                {
                    self.genextend(&mut gx, ne + xc, rigidx, xlbx, xubx, ecount, data);
                }
            }
        }
    }

    fn init_extend(
        &mut self,
        g: &mut Graph,
        xlb: usize,
        xub: usize,
        data: &mut Vec<LevelData>,
    ) -> Option<(usize, usize, usize, usize, usize, usize, usize, Vec<usize>)> {
        let n = g.n();
        let deg = (0..n).map(|i| g.deg(i)).collect_vec();
        let mut xlb = xlb;
        self.nodes[n] += 1;
        let nx = n + 1;
        let dmax = deg[n - 1];
        let dcrit: Option<usize> = (self.mindeg + n).checked_sub(self.maxn);
        let (d, dlow) = d_dlow(n, deg, dmax, dcrit);
        if xlb == dmax && d.count_ones() as usize + dmax > n {
            xlb += 1;
        }
        if nx == self.maxn && xlb < self.mindeg {
            xlb = self.mindeg;
        }
        if xlb > xub {
            return None;
        }
        let imin = data[n].xstart[xlb];
        let imax = data[n].xstart[xub + 1];
        let xorb = data[n].xorb.clone();
        Some((n, nx, dmax, d, dlow, imin, imax, xorb))
    }
}

fn d_dlow(n: usize, deg: Vec<usize>, dmax: usize, dcrit: Option<usize>) -> (usize, usize) {
    let mut d: usize = 0;
    let mut dlow = 0;
    for i in 0..n {
        if deg[i] == dmax {
            d |= xbit(i);
        }
        if Some(true) == dcrit.map(|dcrit| dcrit == deg[i]) {
            dlow |= xbit(i);
        }
    }
    (d, dlow)
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
    let n = g.n();
    if n == 1 {
        *code = 1;
        return;
    }
    *code = 0;
    let mut lact = active.clone();
    while !partition.is_discrete() && !lact.is_zero() {
        let split1 = lact.takebit();
        let split2 = partition.cell_end(split1);
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
    use crate::nauty::{partition_nest::PartitionNest, u32_to_bitvec};
    use bitvec::{order::Msb0, vec::BitVec};
    use test_case::test_case;

    #[test_case(Graph::from_edges(7, &[(0, 5),(0, 6),(1, 6)]), &[1,2,3,4,5,0,6], &[1,1,1,1,0,1,0], 0, 2, u32_to_bitvec(2214592512, 7), false, 0, &[2,3,4,1,5,0,6], &[1,1,0,1,0,1,0], 3, 0; "F?AE?")]
    #[test_case(Graph::from_edges(6, &[(0, 5)]), &[1,2,3,4,5,0], &[1,1,1,0,1,0], 0, 2, u32_to_bitvec(2281701376, 6), false, 0, &[1,2,3,4,5,0], &[1,1,1,0,1,0], 2, 0; "E?A?")]
    #[test_case(Graph::from_edges(7, &[(0, 6)]), &[1,2,3,4,5,0,6], &[1,1,1,1,0,1,0], 0, 2, u32_to_bitvec(2214592512, 7), true, 0, &[1,2,3,4,5,0,6], &[1,1,1,1,0,1,0], 2, 0; "F??C?")]
    #[test_case(Graph::no_edge(3), &[2,1,0], &[1,1,0], 0, 1, u32_to_bitvec(2147483648, 3), false, 1, &[2,1,0], &[1,1,0], 1, 0)]
    #[test_case(Graph::no_edge(2), &[1,0], &[1, 0], 0, 1, u32_to_bitvec(2147483648, 2), false, 32767, &[1, 0], &[1,0], 1, 0)]
    #[allow(clippy::too_many_arguments)]
    fn test_refinex(
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
