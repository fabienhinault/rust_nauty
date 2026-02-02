use super::g6string::BIAS6;

/// G6Char represents a char in a graph6 representation of a graph.
pub struct G6Char {
    value_so_far: u8,
    push_left_count: u8,
}

impl G6Char {
    pub fn new() -> Self {
        Self {
            value_so_far: 0,
            push_left_count: 6,
        }
    }

    pub fn push(&mut self, x: bool) {
        self.value_so_far <<= 1;
        if x {
            self.value_so_far |= 1;
        }
        self.push_left_count -= 1;
    }

    pub fn is_complete(&self) -> bool {
        self.push_left_count == 0
    }

    pub fn is_new(&self) -> bool {
        self.push_left_count == 6
    }

    pub fn to_char(&self) -> u8 {
        BIAS6 + (self.value_so_far << self.push_left_count)
    }
}
