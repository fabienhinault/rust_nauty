use crate::nauty::graph;
use clap::Parser;
use cli::GengCli;
use env::Env;
use geng::Extender;
use level_data::LevelData;
use nauty::{WORDSIZE, graph::Graph};
use once_cell::sync::Lazy;
use std::sync::LazyLock;

pub mod cli;
mod env;
mod geng;
mod graph6;
pub mod gtools;
mod level_data;
mod naugraph;
mod nautil;
pub mod nauty;

const MAXN: usize = WORDSIZE;
static PARSED_ARGS: Lazy<GengCli> = Lazy::new(GengCli::parse);
static ARG_MAXN: LazyLock<usize> = LazyLock::new(|| PARSED_ARGS.maxn);
static MAXE: LazyLock<usize> = LazyLock::new(|| {
    PARSED_ARGS
        .maxe
        .unwrap_or((*ARG_MAXN * *ARG_MAXN - *ARG_MAXN) / 2)
});
static SPARSE: Lazy<bool> = Lazy::new(|| {
    let mut sparse = PARSED_ARGS.bipartite
        || PARSED_ARGS.squarefree
        || PARSED_ARGS.trianglefree
        || PARSED_ARGS.savemem;
    if *ARG_MAXN > 28 || *ARG_MAXN + 4 > 8 * 64usize {
        sparse = true;
    }

    if *ARG_MAXN == *MAXE + 1 && (PARSED_ARGS.connec1 || PARSED_ARGS.connec2) {
        sparse = true;
    }
    sparse
});
static MAX_DEG: Lazy<usize> = Lazy::new(|| {
    let mut maxdeg = PARSED_ARGS.maxdeg.unwrap_or(MAXN);
    if maxdeg >= *MAXE {
        maxdeg = *ARG_MAXN - 1;
    }
    maxdeg
});
static DATA: LazyLock<Vec<LevelData>> =
    LazyLock::new(|| LevelData::make(*SPARSE, *ARG_MAXN, *MAX_DEG));
fn main() {
    let mut args = PARSED_ARGS.clone();
    let maxn = args.maxn;
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

    let mut _g = graph::Graph::default();
    let mut _deg: [usize; 1] = [0];

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
        println!("{}", Graph::one().to_graph6());
    } else {
        let env = Env {
            maxn: 2,
            maxe: 1,
            mine: 0,
            maxdeg: 1,
        };
        let mut sparse = args.bipartite || args.squarefree || args.trianglefree || args.savemem;
        if maxn > 28 || maxn + 4 > 8 * 64usize {
            args.savemem = true;
            sparse = true;
        }
        if maxn == maxe + 1 && connec > 0 {
            args.bipartite = true;
            args.squarefree = true;
            sparse = true;
        }
        let mut data: Vec<LevelData> = LevelData::make(sparse, maxn, maxdeg);
        data[1].xbnds(&env, 1, 0, 0);
        let xlb = data[1].xlb;
        let xub = data[1].xub;
        let mut extender = Extender {
            maxn,
            mindeg,
            nodes: vec![],
            connec,
            canonise: args.canonise,
            args,
        };
        let mut deg = [0];
        let mut g = Graph::no_edge(0);
        /* counts by number of edges */
        let mut ecount: Vec<usize> = vec![];
        extender.genextend(&mut g, 1, &mut deg, 0, true, xlb, xub, &mut ecount, DATA);
    }
}
