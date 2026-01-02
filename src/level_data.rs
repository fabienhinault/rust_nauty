use std::cmp::{max, min};

use crate::{env::Env, nauty::SetWord};

const WORD_SIZE: usize = size_of::<SetWord>();
const MAXN: usize = WORD_SIZE;

#[derive(Default)]
pub struct LevelData {
    edge_count: usize,
    dmax: usize,
    xlb: usize,
    xub: usize,
    lo: SetWord,
    hi: SetWord,
    xstart: [SetWord; MAXN + 1],
    xset: Vec<SetWord>,
    xcard: Vec<SetWord>,
    xinv: Vec<SetWord>,
    xorb: Vec<SetWord>,
    xx: Vec<SetWord>,
    xlim: SetWord,
}

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
