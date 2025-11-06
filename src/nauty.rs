use crate::graph6::{self, encode_graph_size, g6_len, BIAS6};

pub struct Graph {
    pub set_count: usize,
    pub sets: Vec<Set>,
}

pub struct Set {
    pub word_count: usize,
    pub words: Vec<SetWord>,
}

pub type SetWord = u64;

impl Graph {
    pub fn to_graph6(&self) -> String {
        let mut gcode = Vec::with_capacity(g6_len(self.sets.len()));
        gcode.append(&mut encode_graph_size(self.sets.len()));
        let mut k = 6;
        let mut x: u8 = 0;
        for (j, row) in self.sets.iter().enumerate() {
            for i in 0..j {
                x <<= 1;
                if row.is_element(i) {
                    x |= 1;
                    k -= 1;
                    if k == 0 {
                        gcode.push(BIAS6 + x);
                    }
                }
            }
        }
        if k != 6 {
            gcode.push(BIAS6 + x);
        }
        gcode.push('\n' as u8);
        String::from_utf8(gcode).expect("String::from_utf8")
    }
}

impl Set {
    pub fn is_element(&self, pos: usize) -> bool {
        self.words[Self::set_wd(pos)] & graph6::BIT[Self::set_bt(pos)] != 0
    }

    pub fn set_wd(pos: usize) -> usize {
        pos >> 6
    }

    pub fn set_bt(pos: usize) -> usize {
        pos & 0x3F
    }
}