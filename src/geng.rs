use crate::{
    graph::Graph,
    level_data::LevelData,
    nauty::{NautyCounter, graph},
};

#[inline(always)]
fn xbit(i: usize) -> usize {
    1 << i
}

fn xnextbit(i: usize) -> usize {
    i.leading_zeros() as usize
}
pub struct Extender {
    pub maxn: usize,
    pub mindeg: usize,
    pub connec: u8,
    pub nodes: Vec<NautyCounter>,
    pub data: Vec<LevelData>,
    pub canonise: bool,
}

impl Extender {
    /* extend from n to n+1 -- version for general graphs */
    // 2314
    pub fn genextend(
        &mut self,
        g: &mut graph::Graph,
        n: usize,
        deg: &[usize],
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
                if accept2(
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
                if accept1(
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
        }
    }
}

fn accept1(
    g: &mut graph::Graph,
    n: usize,
    x: usize,
    gx: &mut graph::Graph,
    deg: &[usize],
    d: bool,
) -> bool {
    todo!()
}

fn accept2(
    g: &mut graph::Graph,
    n: usize,
    x: usize,
    gx: &mut graph::Graph,
    deg: &[usize],
    d: bool,
) -> bool {
    todo!()
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
