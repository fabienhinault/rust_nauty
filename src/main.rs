use clap::Parser;
use cli::GengCli;
use env::Env;
use geng::Extender;
use level_data::{LevelData, xbnds};
use nauty::{Graph, Set, SetWord, WORDSIZE};

pub mod cli;
mod env;
mod geng;
mod graph6;
pub mod gtools;
mod level_data;
pub mod nauty;

const MAXN: usize = WORDSIZE;
fn main() {
    let parsed_args = GengCli::parse();
    let mut args = parsed_args.clone();
    let maxn = args.maxn as usize;
    let maxe = args.maxe.unwrap_or((maxn * maxn - maxn) / 2);
    let mut maxdeg = args.maxdeg.unwrap_or(MAXN);
    let mut mindeg = args.mindeg.unwrap_or(0);

    if args.connec1 && mindeg < 1 && maxn > 1 {
        mindeg = 1;
    }
    if args.connec2 && mindeg < 2 && maxn > 2 {
        mindeg = 2;
    }
    if maxdeg >= maxe {
        maxdeg = maxn - 1;
    }
    let mut ecount = vec![0; 1 + maxn * (maxn - 1)];

    let mut g = Graph(vec![Set { words: [0] }]);
    let mut deg: [usize; 1] = [0];

    let connec: u8;
    if args.connec2 {
        connec = 2;
    } else if args.connec1 {
        connec = 1;
    } else {
        connec = 0;
    }
    if maxn == 1 {
        ecount[0] += 1;
        println!("{}", g.to_graph6());
    } else {
        let env = Env {
            maxn: 2,
            maxe: 1,
            mine: 0,
            maxdeg: 1,
        };
        let mut sparse = args.bipartite || args.squarefree || args.trianglefree || args.savemem;
        if maxn > 28 || maxn + 4 > 8 * SetWord::BITS as usize {
            args.savemem = true;
            sparse = true;
        }
        if maxn == maxe + 1 && connec > 0 {
            args.bipartite = true;
            args.squarefree = true;
            sparse = true;
        }
        let mut data: Vec<LevelData> = LevelData::make(sparse, maxn, maxdeg);
        xbnds(&env, 1, 0, 0, &mut data[1]);
        let xlb = data[1].xlb;
        let xub = data[1].xub;
        let mut extender = Extender {
            maxn,
            mindeg,
            nodes: vec![],
            data,
        };
        extender.genextend(&mut g, 1, &deg, 0, true, xlb, xub);
    }
}
