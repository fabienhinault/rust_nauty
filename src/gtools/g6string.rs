use std::fmt::Display;

use crate::nauty::Graph;

use super::g6char::G6Char;

pub struct G6String {
    so_far: Vec<u8>,
    current_char: G6Char,
}

const SMALL_N: usize = 62;
const SMALLISH_N: usize = 258047;
pub const BIAS6: u8 = 63;
const MAXBYTE: u8 = 126;
const C6MASK: usize = 63;

impl G6String {
    pub fn from<const M: usize>(g: &Graph<M>) -> Self {
        let mut g6_string = Self::new(g.n());
        for (i_row, row) in g.0.iter().enumerate() {
            for i_other_vertex in 0..i_row {
                g6_string.push(row.is_element(i_other_vertex));
            }
        }
        g6_string
    }

    pub fn new(vertex_number: usize) -> Self {
        let mut so_far = encode_graph_size(vertex_number);
        so_far.reserve(g6_body_len(vertex_number));
        Self {
            so_far,
            current_char: G6Char::new(),
        }
    }

    pub fn push(&mut self, is_element_i_j: bool) {
        self.current_char.push(is_element_i_j);
        if self.current_char.is_complete() {
            self.so_far.push(self.current_char.to_char());
            self.current_char = G6Char::new();
        }
    }

    pub fn finish(&mut self) {
        if !self.current_char.is_new() {
            self.so_far.push(self.current_char.to_char());
        }
        self.so_far.push(b'\n');
    }
}

impl Display for G6String {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            String::from_utf8(self.so_far.clone()).expect("String::from_utf8")
        )
    }
}

// function N(n) in https://users.cecs.anu.edu.au/~bdm/data/formats.txt
pub fn encode_graph_size(vertex_number: usize) -> Vec<u8> {
    if vertex_number <= SMALL_N {
        vec![BIAS6 + vertex_number as u8]
    } else if vertex_number <= SMALLISH_N {
        vec![
            MAXBYTE,
            BIAS6 + ((vertex_number >> 12) as u8),
            BIAS6 + ((vertex_number >> 6) & C6MASK) as u8,
            BIAS6 + (vertex_number & C6MASK) as u8,
        ]
    } else {
        vec![
            MAXBYTE,
            MAXBYTE,
            BIAS6 + ((vertex_number >> 30) as u8),
            BIAS6 + ((vertex_number >> 24) & C6MASK) as u8,
            BIAS6 + ((vertex_number >> 18) & C6MASK) as u8,
            BIAS6 + ((vertex_number >> 12) & C6MASK) as u8,
            BIAS6 + ((vertex_number >> 6) & C6MASK) as u8,
            BIAS6 + (vertex_number & C6MASK) as u8,
        ]
    }
}

fn g6_body_len(n: usize) -> usize {
    (n / 12) * (n - 1) + ((n % 12) * (n - 1) + 11) / 12
}
