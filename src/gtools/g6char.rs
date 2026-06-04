use crate::gtools::g6error::G6Error;

use super::g6string::BIAS6;

const TOPBIT6: u8 = 32;
/// G6Char represents a char in a graph6 representation of a graph.
pub struct G6Char {
    value_so_far: u8,
    push_left_count: u8,
}

impl G6Char {
    const CAPACITY: u8 = 6;
    pub fn new() -> Self {
        Self {
            value_so_far: 0,
            push_left_count: Self::CAPACITY,
        }
    }

    pub fn push_back(&mut self, x: bool) {
        self.value_so_far <<= 1;
        if x {
            self.value_so_far |= 1;
        }
        self.push_left_count -= 1;
    }

    pub fn pop_front(&mut self) -> bool {
        let top_bit = self.value_so_far & TOPBIT6 != 0;
        self.value_so_far <<= 1;
        self.push_left_count += 1;
        top_bit
    }

    pub fn len(&self) -> u8 {
        Self::CAPACITY - self.push_left_count
    }

    pub fn is_full(&self) -> bool {
        self.len() == Self::CAPACITY
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn to_char(&self) -> u8 {
        BIAS6 + (self.value_so_far << self.push_left_count)
    }
}

impl TryFrom<u8> for G6Char {
    type Error = G6Error;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        Ok(Self {
            value_so_far: value.checked_sub(BIAS6).ok_or(G6Error())?,
            push_left_count: 0,
        })
    }
}
