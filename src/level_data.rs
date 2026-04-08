use crate::env::Env;
//use bitvec::{order::Msb0, vec::BitVec};
use std::cmp::min;

const WORD_SIZE: usize = size_of::<usize>();
const MAXN: usize = WORD_SIZE;

#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub struct SetCard {
    card: usize, /* cardinalities of all x-sets */
    //set: BitVec<usize, Msb0>, /* array of all x-sets in card order */
    set: usize,
}
#[derive(Default)]
pub struct LevelData {
    edge_count: usize, /* values used for xlb,xub calculation */
    dmax: usize,
    /* saved bounds on extension degree */
    // extension degree lower bound
    pub xlb: usize,
    // extension degree upper bound
    pub xub: usize,
    lo: usize, /* work purposes for orbit calculation */
    hi: usize,
    pub xstart: [usize; MAXN + 1], /* index into xset[] for each cardinality */
    // xset: Vec<SetWord>,
    // xcard: Vec<SetWord>,
    x_set_card: Vec<SetCard>, /* array of all x-sets in card order, cardinalities of all x-sets */
    xinv: Vec<Option<usize>>, /* map from x-set to index in xset */
    pub xorb: Vec<usize>,     /* min orbit representative */
    xx: Vec<usize>,           /* (-b, -t, -s, -m) candidate x-sets */
    /*   note: can be the same as xcard */
    xlim: usize, /* number of x-sets in xx[] */
}

fn arith(a: usize, b: usize, c: usize) -> usize {
    (a / c) * b + ((a % c) * b) / c
}

impl LevelData {
    pub fn xcard(&self, i: usize) -> usize {
        self.x_set_card[i].card
    }

    // pub fn xset(&self, i: usize) -> BitVec<usize, Msb0> {
    //     self.x_set_card[i].set.clone()
    // }

    pub fn xset(&self, i: usize) -> usize {
        self.x_set_card[i].set
    }

    //static void makeleveldata(boolean restricted)
    /* make the level data for each level */
    pub fn make(restricted: bool, maxn: usize, maxdeg: usize) -> Vec<Self> {
        let mut data: Vec<LevelData> = Vec::with_capacity(maxn);
        data.push(LevelData::default());
        for n in 1..maxn {
            let nn = min(n, maxdeg);
            let mut ncj: usize = 1;
            let mut nxsets: usize = 1;
            for j in 1..nn {
                ncj = arith(ncj, n - j + 1, j);
                nxsets += ncj;
            }
            let mut d = Self::default();
            d.xlb = usize::MAX;

            if restricted {
                d.xorb = Vec::with_capacity(nxsets);
                d.xx = Vec::with_capacity(nxsets);
                continue;
            }
            let tttn: usize = 1 << n;
            d.x_set_card = Vec::with_capacity(nxsets);
            d.xinv = vec![None; tttn];
            d.xorb = Vec::with_capacity(nxsets);

            let ilast = if n == WORD_SIZE {
                usize::MAX
            } else {
                (1 << n) - 1
            };
            for i in 0..ilast {
                let h: usize = i.count_ones() as usize;
                if h <= maxdeg {
                    d.x_set_card.push(SetCard { card: h, set: i });
                }
                if i == ilast {
                    break;
                }
            }

            if d.x_set_card.len() != nxsets {
                panic!();
            }

            d.x_set_card.sort();

            for i in 0..nxsets {
                d.xinv[d.x_set_card[i].set] = Some(i);
            }

            d.xstart[0] = 0;
            for i in 1..nxsets {
                if d.x_set_card[i].card > d.x_set_card[i - 1].card {
                    d.xstart[d.xcard(i) as usize] = i;
                }
            }
            d.xstart[(d.xcard(nxsets - 1) + 1) as usize] = nxsets;
            data.push(d);
        }

        // Initialize xstart arrays for each level
        for level in &mut data {
            level.xstart = [0; MAXN + 1];
            level.x_set_card = Vec::new();
            level.xinv = Vec::new();
            level.xorb = Vec::new();
            level.xx = Vec::new();
        }

        data
    }

    /* find bounds on extension degree;  store answer in data[*].*  */
    pub fn xbnds(&mut self, env: &Env, n: usize, edge_count: usize, dmax: usize) {
        let mut xlb = if n == 1 {
            0
        } else if dmax > (2 * edge_count + n - 2) / (n - 1) {
            dmax
        } else {
            (2 * edge_count + n - 2) / (n - 1)
        };
        let mut xub = if n < MAXN { n } else { MAXN };
        let mut d;
        let mut m;
        for xc in (xlb..=xub).rev() {
            d = xc;
            m = edge_count + d;
            for nn in n + 1..env.maxn {
                if d < (2 * m + nn - 2) / (nn - 1) {
                    d = (2 * m + nn - 2) / (nn - 1);
                }
                m += d;
            }
            if d > env.maxdeg || m > env.maxe {
                xub = xc - 1;
            } else {
                break;
            }
        }

        if edge_count + xlb < env.mine {
            for xc in xlb..=xub {
                m = edge_count + xc;
                for _nn in n + 1..env.maxn {
                    m += min(n, env.maxdeg);
                }
                if m < env.mine {
                    xlb = xc + 1;
                } else {
                    break;
                }
            }
        }
        self.edge_count = edge_count;
        self.dmax = dmax;
        self.xlb = xlb;
        self.xub = xub;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_arith() {
        let a: usize = usize::MAX / 2;
        let b: usize = 4;
        let c: usize = 4;
        // let d = a * b / c; // error: this operation will overflow
        let d = arith(a, b, c);
        assert_eq!(d, a);
    }
}
