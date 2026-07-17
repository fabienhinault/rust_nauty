use super::*;
use test_case::test_case;

pub mod graphs;

#[test]
fn test_bitvec() {
    let bv = bitvec![usize, Msb0; 1, 0, 0];
    assert!(bv[0]);
    assert!(!bv[1]);
    assert!(!bv[2]);
}

// #[test]
// fn test_bit() {
//     // for i in 0.. 10000 {
//     //     let j =
//     // }
//     assert_eq!(64, BIT.len());
//     for (i_bit, bit) in BIT.into_iter().enumerate() {
//         assert_eq!(1 << (63 - i_bit), bit)
//     }
// }

// example from https://users.cecs.anu.edu.au/~bdm/data/formats.txt, line 73
#[test]
fn test_to_graph6() {
    let g = create_example();
    let g6 = g.to_graph6();
    assert_eq!(g6.bytes().collect::<Vec<_>>(), [68, 81, 99, 10]);
    assert_eq!(g6, "DQc\n");
}

// example from https://users.cecs.anu.edu.au/~bdm/data/formats.txt, line 73
#[test]
fn test_from_graph6() {
    assert_eq!(
        create_example(),
        Graph::from_graph6("DQc".to_owned()).unwrap()
    );
    assert_eq!(
        create_example(),
        Graph::from_graph6("DQc\n".to_owned()).unwrap()
    );
}

#[test]
fn test_example_is_connected() {
    assert!(create_example().isconnected());
}

#[test]
fn test_disconnected() {
    assert!(!create_disconnected().isconnected());
}

#[test]
fn test_disconnected_isbiconnected() {
    assert!(!create_disconnected().isbiconnected());
}

#[test]
fn test_example_isbiconnected() {
    assert!(!create_example().isbiconnected());
}

#[test_case(graphs::cycle(6))]
#[test_case(create_diamond())]
#[test_case(create_complete(4))]
#[test_case(create_g4g_biconnected())]
fn test_is_biconnected(biconnected_graph: Graph) {
    assert!(biconnected_graph.isbiconnected());
}

#[test_case(graphs::path(4))]
#[test_case(create_g4g_not_biconnected())]
#[test_case(create_D_7dC())]
fn test_is_not_biconnected(not_biconnected_graph: Graph) {
    // println!("{}", not_biconnected_graph.to_graph6());
    assert!(!not_biconnected_graph.isbiconnected());
}

#[test]
fn test_diamond_isbiconnected() {
    assert!(create_diamond().isbiconnected());
}

#[test]
fn test_tetraedron_isbiconnected() {
    assert!(create_complete(4).isbiconnected());
}

#[test]
fn test_g4g_not_isbiconnected() {
    assert!(!create_g4g_not_biconnected().isbiconnected());
}

#[test]
fn test_without_loop_not_isbiconnected() {
    assert!(!create_without_loop().isbiconnected());
}

#[test]
fn test_create_from_u32() {
    let actual = Graph::from_u32(&[1610612736, 2684354560, 3221225472]);
    let expected = Graph(vec![
        //                   0  1  2
        bitvec![usize, Msb0; 0, 1, 1],
        bitvec![usize, Msb0; 1, 0, 1],
        bitvec![usize, Msb0; 1, 1, 0],
    ]);
    assert_eq!(actual, expected);
}

#[test]
fn test_complete_bipartite() {
    println!("{}", graphs::complete_bipartite(2, 3).to_matrix());
    assert_eq!(
        graphs::complete_bipartite(2, 3).to_matrix(),
        "0 0 1 1 1\n0 0 1 1 1\n1 1 0 0 0\n1 1 0 0 0\n1 1 0 0 0"
    );
}

#[test]
fn test_complete_tripartite() {
    println!("{}", graphs::complete_tripartite(2, 2, 2).to_matrix());
    assert_eq!(
        graphs::complete_tripartite(2, 2, 2).to_matrix(),
        "0 0 1 1 1 1\n0 0 1 1 1 1\n1 1 0 0 1 1\n1 1 0 0 1 1\n1 1 1 1 0 0\n1 1 1 1 0 0"
    );
}

#[test]
fn test_complete_k_partite() {
    println!("{}", graphs::complete_k_partite(&[2, 2, 2, 2]).to_matrix());
    assert_eq!(
        graphs::complete_k_partite(&[2, 2, 2, 2]).to_matrix(),
        "0 0 1 1 1 1 1 1\n0 0 1 1 1 1 1 1\n1 1 0 0 1 1 1 1\n1 1 0 0 1 1 1 1\n1 1 1 1 0 0 1 1\n1 1 1 1 0 0 1 1\n1 1 1 1 1 1 0 0\n1 1 1 1 1 1 0 0"
    );
}

#[test]
fn test_wheel() {
    println!("{}", graphs::wheel(6).to_matrix());
    assert_eq!(
        graphs::wheel(6).to_matrix(),
        "0 1 1 1 1 1\n1 0 1 0 0 1\n1 1 0 1 0 0\n1 0 1 0 1 0\n1 0 0 1 0 1\n1 1 0 0 1 0"
    );
}

#[test]
fn test_join() {
    let g = Graph::join(&Graph::one(), &graphs::path(3));
    println!("{}", g.to_matrix());
    assert_eq!(g, create_diamond());
}

// example from https://users.cecs.anu.edu.au/~bdm/data/formats.txt, line 73
//
//  2---0---4---3---1
//
fn create_example() -> Graph {
    Graph(vec![
        //                   0  1  2  3  4
        bitvec![usize, Msb0; 0, 0, 1, 0, 1],
        bitvec![usize, Msb0; 0, 0, 0, 1, 0],
        bitvec![usize, Msb0; 1, 0, 0, 0, 0],
        bitvec![usize, Msb0; 0, 1, 0, 0, 1],
        bitvec![usize, Msb0; 1, 0, 0, 1, 0],
    ])
}

//  2---0   4---3---1
fn create_disconnected() -> Graph {
    Graph(vec![
        //                   0  1  2  3  4
        bitvec![usize, Msb0; 0, 0, 1, 0, 0],
        bitvec![usize, Msb0; 0, 0, 0, 1, 0],
        bitvec![usize, Msb0; 1, 0, 0, 0, 0],
        bitvec![usize, Msb0; 0, 1, 0, 0, 1],
        bitvec![usize, Msb0; 0, 0, 0, 1, 0],
    ])
}

pub fn create_f(n: usize, f: fn(usize, usize, usize) -> bool) -> Graph {
    Graph((0..n).map(|i| bitvec_from_f(n, i, f)).collect())
}

pub fn create_star() -> Graph {
    create_f(5, |i_current_vertex, i_other_vertex, _n| {
        i_other_vertex == (i_current_vertex + 3).rem_euclid(5)
            || i_other_vertex == (i_current_vertex + 2).rem_euclid(5)
    })
}

//   1
//  / \
// 0---2
//  \ /
//   3
pub fn create_diamond() -> Graph {
    Graph(vec![
        //                   0  1  2  3
        bitvec![usize, Msb0; 0, 1, 1, 1],
        bitvec![usize, Msb0; 1, 0, 1, 0],
        bitvec![usize, Msb0; 1, 1, 0, 1],
        bitvec![usize, Msb0; 1, 0, 1, 0],
    ])
}

pub fn create_complete(n: usize) -> Graph {
    Graph((0..n).map(|i| create_complete_bitvec(n, i)).collect())
}

pub fn create_zero(n: usize) -> Graph {
    Graph::no_edge(n)
}

fn create_complete_bitvec(n: usize, i: usize) -> Set {
    let mut result: BitVec<usize, Msb0> = bitvec![usize, Msb0; 1; n];
    result.set(i, false);
    result
}

// https://www.geeksforgeeks.org/dsa/biconnectivity-in-a-graph/
// 1-0--3
// |/   |
// 2    4
fn create_g4g_not_biconnected() -> Graph {
    Graph(vec![
        //                   0  1  2  3  4
        bitvec![usize, Msb0; 0, 1, 1, 1, 0],
        bitvec![usize, Msb0; 1, 0, 1, 0, 0],
        bitvec![usize, Msb0; 1, 1, 0, 0, 0],
        bitvec![usize, Msb0; 1, 0, 0, 0, 1],
        bitvec![usize, Msb0; 0, 0, 0, 1, 0],
    ])
}

//  1-0--3
//  |/   |
//  2----4
fn create_g4g_biconnected() -> Graph {
    Graph(vec![
        //                   0  1  2  3  4
        bitvec![usize, Msb0; 0, 1, 1, 1, 0],
        bitvec![usize, Msb0; 1, 0, 1, 0, 0],
        bitvec![usize, Msb0; 1, 1, 0, 0, 1],
        bitvec![usize, Msb0; 1, 0, 0, 0, 1],
        bitvec![usize, Msb0; 0, 0, 1, 1, 0],
    ])
}

//  ,----.
//  1-0--3
//  |/   |
//  2    4
//
// D}C
fn create_D_7dC() -> Graph {
    Graph(vec![
        //                   0  1  2  3  4
        bitvec![usize, Msb0; 0, 1, 1, 1, 0],
        bitvec![usize, Msb0; 1, 0, 1, 1, 0],
        bitvec![usize, Msb0; 1, 1, 0, 0, 0],
        bitvec![usize, Msb0; 1, 1, 0, 0, 1],
        bitvec![usize, Msb0; 0, 0, 0, 1, 0],
    ])
}

fn create_without_loop() -> Graph {
    Graph(vec![
        //                   0  1  2  3
        bitvec![usize, Msb0; 0, 1, 1, 1],
        bitvec![usize, Msb0; 1, 0, 0, 0],
        bitvec![usize, Msb0; 1, 0, 0, 0],
        bitvec![usize, Msb0; 1, 0, 0, 0],
    ])
}
