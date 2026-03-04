use crate::{Graph, nauty::SetWord};

pub struct Extender {
    maxn: usize,
    mindeg: usize,
    dlow: usize,
}

impl Extender {
    fn genextend(
        &self,
        g: &mut Graph<1>,
        n: usize,
        deg: &[usize],
        ne: usize,
        rigid: bool,
        xlb: usize,
        xub: usize,
    ) {
        let nx = n + 1;
        let dmax = deg[n - 1];
        let dcrit = self.mindeg - self.maxn + n;
        let mut d = 0;
        let mut dlow = 0;
        for i in 0..n {
            if deg[i] == dmax {
                d |= xbit(i);
            }
            if deg[i] == dcrit {
                dlow |= xbit(i);
            }
        }
    }
}

#[inline(always)]
fn xbit(i: usize) -> SetWord {
    1 << i
}
