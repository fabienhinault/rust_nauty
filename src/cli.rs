use crate::nauty::WORDSIZE;
use clap::Parser;

#[derive(Parser, Clone)]
#[command(version)]
pub struct GengCli {
    #[arg(short = 'n', default_value_t = false)]
    pub nautyformat: bool,
    #[arg(short = 'u', default_value_t = false)]
    pub nooutput: bool,
    #[arg(short = 'g', default_value_t = true)]
    pub graph6: bool,
    #[arg(short = 's', default_value_t = false)]
    pub sparse6: bool,
    #[arg(short = 't', default_value_t = false)]
    pub trianglefree: bool,
    #[arg(short = 'f', default_value_t = false)]
    pub squarefree: bool,
    #[arg(short = 'b', default_value_t = false)]
    pub bipartite: bool,
    #[arg(short = 'v', default_value_t = false)]
    verbose: bool,
    #[arg(short = 'l', default_value_t = false)]
    canonise: bool,
    #[arg(short = 'y', default_value_t = false)]
    yformat: bool,
    // diff with original geng: -h is for help
    #[arg(long, default_value_t = false)]
    header: bool,
    #[arg(short = 'm', default_value_t = false)]
    pub savemem: bool,
    #[arg(short = 'c', default_value_t = false)]
    pub connec1: bool,
    #[arg(short = 'C', default_value_t = false)]
    pub connec2: bool,
    #[arg(short = 'q', default_value_t = false)]
    quiet: bool,
    #[arg(short = '$', default_value_t = false)]
    secret: bool,
    #[arg(short = 'S', default_value_t = false)]
    safe: bool,
    #[arg(short = 'd')]
    pub mindeg: Option<usize>,
    #[arg(short = 'D')]
    pub maxdeg: Option<usize>,
    #[arg(short = 'x')]
    multiplicity: Option<i32>,
    // split level increment
    #[arg(short = 'X')]
    splitlevinc: Option<i32>,
    // diff with original geng
    #[arg(long)]
    pub maxe: Option<usize>,
    #[arg(value_parser = clap::value_parser!(u64).range(1..=(WORDSIZE as u64)))]
    pub maxn: u64,
}
