mod env;
mod graph6;
mod level_data;
mod nauty;

use env::Env;
use level_data::{LevelData, xbnds};
use nauty::{Graph, Set};

fn main() {
    let env = Env {
        maxn: 2,
        maxe: 1,
        mine: 0,
        maxdeg: 1,
    };
    let mut data: LevelData = LevelData::default();
    xbnds(&env, 1, 0, 0, &mut data);
    let g = Graph(vec![Set { words: [0] }]);
    println!("{}", g.to_graph6());
}
