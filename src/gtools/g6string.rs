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
    pub fn from(g: &Graph) -> Self {
        g.into()
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
        self.current_char.push_back(is_element_i_j);
        if self.current_char.is_full() {
            self.so_far.push(self.current_char.to_char());
            self.current_char = G6Char::new();
        }
    }

    fn finish(&mut self) {
        if !self.current_char.is_empty() {
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

impl From<&Graph> for G6String {
    fn from(g: &Graph) -> Self {
        let mut g6_string = Self::new(g.n());
        for (i_row, row) in g.0.iter().enumerate() {
            for i_other_vertex in 0..i_row {
                g6_string.push(row[i_other_vertex]);
            }
        }
        g6_string.finish();
        g6_string
    }
}

impl From<&str> for G6String {
    fn from(value: &str) -> Self {
        Self::new(0)
    }
}

// function N(n) in https://users.cecs.anu.edu.au/~bdm/data/formats.txt
// function encodegraphsize in gtools.c
/* Encode the size n in a string starting at **p, and reset **p
to point to the character after the size */
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

fn next(iter: &mut std::slice::Iter<'_, u8>) -> Result<usize, ()> {
    Ok((iter.next().ok_or(())?.checked_sub(BIAS6).ok_or(())?) as usize)
}

// function graphsize in gtools.c
/* Get size of graph out of graph6, digraph6 or sparse6 string. */
pub(crate) fn graph_size(iter: &mut std::slice::Iter<'_, u8>) -> Result<usize, ()> {
    let mut n: usize = next(iter)? as usize;
    if n > SMALL_N {
        n = next(iter)?;
        if n > SMALL_N {
            n = next(iter)?;
            n = (n << 6) | next(iter)?;
            n = (n << 6) | next(iter)?;
            n = (n << 6) | next(iter)?;
            n = (n << 6) | next(iter)?;
            n = (n << 6) | next(iter)?;
        } else {
            n = (n << 6) | next(iter)?;
            n = (n << 6) | next(iter)?;
        }
    }
    Ok(n)
}

fn g6_body_len(n: usize) -> usize {
    (n / 12) * (n - 1) + ((n % 12) * (n - 1)).div_ceil(12)
}
