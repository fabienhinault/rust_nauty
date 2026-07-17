use crate::nauty::Graph;

fn eq_mod(i: usize, j: usize, n: usize) -> bool {
    i.rem_euclid(n) == j.rem_euclid(n)
}

fn are_next(i: usize, j: usize, n: usize) -> bool {
    eq_mod(i, j + 1, n) || eq_mod(i + 1, j, n)
}

//  0---1---2---3---4 ... --n
//  `---------------- ... --'
//
pub fn cycle(n: usize) -> Graph {
    Graph::from_f(n, are_next)
}

//  0---1---2---3---4--- ... ---n
pub fn path(n: usize) -> Graph {
    Graph::from_f(n, |i_cur_vertex, i_other_vertex, n| match i_cur_vertex {
        0 => i_other_vertex == 1,
        i if i == n - 1 => i_other_vertex == n - 2,
        _ => are_next(i_cur_vertex, i_other_vertex, n),
    })
}

pub fn star(n: usize) -> Graph {
    Graph::from_f(n, |i_cur_vertex, i_other_vertex, _n| match i_cur_vertex {
        0 => i_other_vertex != 0,
        _ => i_other_vertex == 0,
    })
}

pub fn complete(n: usize) -> Graph {
    Graph::from_f(n, |i_cur_vertex, i_other_vertex, _n| {
        i_cur_vertex != i_other_vertex
    })
}

pub fn wheel(n: usize) -> Graph {
    Graph::from_f(n, |i_cur_vertex, i_other_vertex, n| match i_cur_vertex {
        0 => i_other_vertex != 0,
        _ => i_other_vertex == 0 || are_next(i_cur_vertex - 1, i_other_vertex - 1, n - 1),
    })
}

pub fn complete_bipartite(p: usize, q: usize) -> Graph {
    Graph::from_closure(p + q, |i_cur_vertex, i_other_vertex, _n| {
        (i_cur_vertex < p) == (i_other_vertex >= p)
    })
}

pub fn complete_tripartite(p: usize, q: usize, r: usize) -> Graph {
    Graph::from_closure(p + q + r, |i_cur_vertex, i_other_vertex, _n| {
        if i_cur_vertex < p {
            i_other_vertex >= p
        } else if i_cur_vertex < p + q {
            i_other_vertex < p || i_other_vertex >= p + q
        } else {
            i_other_vertex < p + q
        }
    })
}

pub fn complete_k_partite(ps: &[usize]) -> Graph {
    let (_, indexes) = ps.iter().fold((0, vec![]), |(mut n, mut v), elt| {
        n += elt;
        v.push(n);
        (n, v)
    });
    Graph::from_closure(ps.iter().sum(), |i_cur_vertex, i_other_vertex, _n| {
        let before = indexes.iter().rfind(|&&i| i <= i_cur_vertex);
        let after = indexes.iter().find(|&&i| i > i_cur_vertex);
        match (before, after) {
            (None, None) => false,
            (None, Some(&after)) => after <= i_other_vertex,
            (Some(&before), None) => i_other_vertex < before,
            (Some(&before), Some(&after)) => i_other_vertex < before || after <= i_other_vertex,
        }
    })
}
