use std::cmp::{max, min};

use crate::{env::Env, nauty::SetWord};

const WORD_SIZE: usize = size_of::<SetWord>();
const MAXN: usize = WORD_SIZE;

#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub struct SetCard {
    card: SetWord,
    set: SetWord,
}
#[derive(Default)]
pub struct LevelData {
    edge_count: usize,
    dmax: usize,
    // extension degree lower bound
    pub xlb: usize,
    // extension degree upper bound
    pub xub: usize,
    lo: SetWord,
    hi: SetWord,
    xstart: [SetWord; MAXN + 1],
    // xset: Vec<SetWord>,
    // xcard: Vec<SetWord>,
    x_set_card: Vec<SetCard>,
    xinv: Vec<usize>,
    xorb: Vec<SetWord>,
    xx: Vec<SetWord>,
    xlim: SetWord,
}

fn arith(a: usize, b: usize, c: usize) -> usize {
    (a / c) * b + ((a % c) * b) / c
}

impl LevelData {
    pub fn xcard(&self, i: usize) -> SetWord {
        self.x_set_card[i].card
    }

    pub fn xset(&self, i: usize) -> SetWord {
        self.x_set_card[i].set
    }

    pub fn make(restricted: bool, maxn: usize, maxdeg: usize) -> Vec<Self> {
        let mut data: Vec<LevelData> = Vec::with_capacity(maxn);
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
            d.xinv = Vec::with_capacity(tttn);
            d.xorb = Vec::with_capacity(nxsets);

            let mut j = 0;
            let ilast = if n == WORD_SIZE {
                SetWord::MAX
            } else {
                (1 << n) - 1
            };
            for i in 0..ilast {
                let h: SetWord = i.count_ones().into();
                if h <= maxdeg as SetWord {
                    d.x_set_card[j] = SetCard { card: h, set: i };
                    j += 1;
                }
                if i == ilast {
                    break;
                }
            }

            if j != nxsets {
                panic!();
            }

            d.x_set_card.sort();

            for i in 0..nxsets {
                d.xinv[d.x_set_card[i].set as usize] = i;
            }

            d.xstart[0] = 0;
            for i in 1..nxsets {
                if d.x_set_card[i].card > d.x_set_card[i - 1].card {
                    d.xstart[d.xcard(i) as usize] = i as SetWord;
                }
            }
            d.xstart[(d.xcard(nxsets - 1) + 1) as usize] = nxsets as SetWord;
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
}

/* find bounds on extension degree;  store answer in data[*].*  */
pub fn xbnds(env: &Env, n: usize, edge_count: usize, dmax: usize, data: &mut LevelData) {
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
    data.edge_count = edge_count;
    data.dmax = dmax;
    data.xlb = xlb;
    data.xub = xub;
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
