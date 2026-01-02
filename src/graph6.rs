use crate::nauty::SetWord;

// https://users.cecs.anu.edu.au/~bdm/data/formats.txt

pub const BIT: [SetWord; 64] = [
    0x8000000000000000,
    0x4000000000000000,
    0x2000000000000000,
    0x1000000000000000,
    0x800000000000000,
    0x400000000000000,
    0x200000000000000,
    0x100000000000000,
    0x80000000000000,
    0x40000000000000,
    0x20000000000000,
    0x10000000000000,
    0x8000000000000,
    0x4000000000000,
    0x2000000000000,
    0x1000000000000,
    0x800000000000,
    0x400000000000,
    0x200000000000,
    0x100000000000,
    0x80000000000,
    0x40000000000,
    0x20000000000,
    0x10000000000,
    0x8000000000,
    0x4000000000,
    0x2000000000,
    0x1000000000,
    0x800000000,
    0x400000000,
    0x200000000,
    0x100000000,
    0x80000000,
    0x40000000,
    0x20000000,
    0x10000000,
    0x8000000,
    0x4000000,
    0x2000000,
    0x1000000,
    0x800000,
    0x400000,
    0x200000,
    0x100000,
    0x80000,
    0x40000,
    0x20000,
    0x10000,
    0x8000,
    0x4000,
    0x2000,
    0x1000,
    0x800,
    0x400,
    0x200,
    0x100,
    0x80,
    0x40,
    0x20,
    0x10,
    0x8,
    0x4,
    0x2,
    0x1,
];

const SMALL_N: usize = 62;
const SMALLISH_N: usize = 258047;
pub const BIAS6: u8 = 63;
const MAXBYTE: u8 = 126;
const C6MASK: usize = 63;

pub fn g6_len(n: usize) -> usize {
    size_len(n) + g6_body_len(n)
}

fn size_len(n: usize) -> usize {
    if n <= SMALL_N {
        1
    } else if n <= SMALLISH_N {
        4
    } else {
        8
    }
}

fn g6_body_len(n: usize) -> usize {
    (n / 12) * (n - 1) + ((n % 12) * (n - 1) + 11) / 12
}

// function N(n) in https://users.cecs.anu.edu.au/~bdm/data/formats.txt
pub fn encode_graph_size(n: usize) -> Vec<u8> {
    if n <= SMALL_N {
        vec![BIAS6 + n as u8]
    } else if n <= SMALLISH_N {
        vec![
            MAXBYTE,
            BIAS6 + ((n >> 12) as u8),
            BIAS6 + ((n >> 6) & C6MASK) as u8,
            BIAS6 + (n & C6MASK) as u8,
        ]
    } else {
        vec![
            MAXBYTE,
            MAXBYTE,
            BIAS6 + ((n >> 30) as u8),
            BIAS6 + ((n >> 24) & C6MASK) as u8,
            BIAS6 + ((n >> 18) & C6MASK) as u8,
            BIAS6 + ((n >> 12) & C6MASK) as u8,
            BIAS6 + ((n >> 6) & C6MASK) as u8,
            BIAS6 + (n & C6MASK) as u8,
        ]
    }
}
