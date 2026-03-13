use crate::{Graph, level_data::LevelData, nauty::NautyCounter};

pub struct Extender {
    pub maxn: usize,
    pub mindeg: usize,
    pub connec: u8,
    pub nodes: Vec<NautyCounter>,
    pub data: Vec<LevelData>,
}

impl Extender {
    /* extend from n to n+1 -- version for general graphs */
    pub fn genextend(
        &mut self,
        g: &mut Graph,
        n: usize,
        deg: &[usize],
        ne: usize,
        rigid: bool,
        xlb: usize,
        xub: usize,
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
        let mut gx: Graph = Graph::default();
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
                    gx,
                    deg,
                    xc > dmax + 1 || (xc == dmax + 1 && (x & d) == 0),
                ) && (self.connec == 0
                    || (self.connec == 1 && gx.isconnected())
                    || (self.connec == 2 && gx.isbiconnected()))
                {
                    ecount[ne + xc] += 1;
                    println!(
                        "{}",
                        if cannonise {
                            gcan.to_graph6()
                        } else {
                            gx.to_graph6()
                        }
                    );
                }
            }
        }
    }
}

#[inline(always)]
fn xbit(i: usize) -> SetWord {
    1 << i
}
