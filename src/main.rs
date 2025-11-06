mod nauty;
mod graph6;

use nauty::{Graph, Set};

fn main() {
    let g = Graph{set_count: 1, sets: vec![Set{word_count: 1, words: vec![0]}]};
    println!("{}", g.to_graph6());
}
