use std::{fmt::Debug, mem};

use crate::nauty::NAUTY_INFINITY;

pub struct PartitionNest {
    lab_ptn: Vec<(usize, usize)>,
    pub lab: Vec<usize>,
    pub ptn: Vec<usize>,
}

impl PartitionNest {
    pub fn new(lab: Vec<usize>, ptn: Vec<usize>) -> Self {
        assert_eq!(lab.len(), ptn.len());
        Self {
            lab_ptn: lab.clone().into_iter().zip(ptn.clone()).collect(),
            lab,
            ptn,
        }
    }

    pub fn swap(&mut self, c1: usize, c2: usize) {
        if c1 != c2 {
            let i = c1.min(c2);
            let a = c1.max(c2);
            let (first, second) = self.lab_ptn.split_at_mut(a);
            mem::swap(&mut first[i].0, &mut second[0].0);
        }
    }

    pub fn partition(&self, level: usize) -> Vec<Vec<usize>> {
        self.lab_ptn
            .chunk_by(|(_lab_a, ptn_a), (_lab_b, _ptn_b)| *ptn_a > level)
            .map(|chunk: &[(usize, usize)]| {
                chunk.iter().map(|(lab, _ptn)| *lab).collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    }

    pub fn partition_string(&self, level: usize) -> String {
        itertools::Itertools::intersperse(
            self.partition(level)
                .into_iter()
                .map(|chunk| format!("{chunk:?}"))
                .map(|s| {
                    s.strip_prefix("[")
                        .unwrap()
                        .strip_suffix("]")
                        .unwrap()
                        .to_owned()
                }),
            " | ".to_owned(),
        )
        .collect()
    }

    pub fn lab(&self) -> impl Iterator<Item = usize> {
        self.lab_ptn.iter().map(|(lab, _ptn)| *lab)
    }

    pub fn ptn(&self) -> impl Iterator<Item = usize> {
        self.lab_ptn.iter().map(|(_lab, ptn)| *ptn)
    }

    pub fn max_level(&self) -> Option<usize> {
        self.ptn().filter(|l| l < &NAUTY_INFINITY).max()
    }
}

impl Debug for PartitionNest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PartitionNest")
            .field("lab_ptn", &self.lab_ptn)
            .finish()?;
        if let Some(max_level) = self.max_level() {
            #[allow(clippy::writeln_empty_string)]
            writeln!(f, "")?;
            for l in 0..=max_level {
                writeln!(f, "{l}:  {}", self.partition_string(l))?;
            }
        }
        Ok(())
    }
}

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionNestChunkBy<'a> {
    lab_slice: &'a [usize],
    ptn_slice: &'a [usize],
    level: usize,
}

impl<'a> PartitionNestChunkBy<'a> {
    pub(super) const fn new(lab_slice: &'a [usize], ptn_slice: &'a [usize], level: usize) -> Self {
        PartitionNestChunkBy {
            lab_slice,
            ptn_slice,
            level,
        }
    }
}

impl<'a> Iterator for PartitionNestChunkBy<'a> {
    type Item = &'a [usize];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.lab_slice.is_empty() {
            None
        } else {
            let mut len = 1;
            let mut ptn_iter = self.ptn_slice.iter();
            while let Some(ptn) = ptn_iter.next() {
                if *ptn > self.level {
                    len += 1
                } else {
                    break;
                }
            }
            let (lab_head, lab_tail) = self.lab_slice.split_at(len);
            self.lab_slice = lab_tail;
            let (_ptn_head, ptn_tail) = self.ptn_slice.split_at(len);
            self.ptn_slice = ptn_tail;
            Some(lab_head)
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.lab_slice.is_empty() {
            (0, Some(0))
        } else {
            (1, Some(self.lab_slice.len()))
        }
    }

    #[inline]
    fn last(mut self) -> Option<Self::Item> {
        self.next_back()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::nauty::NAUTY_INFINITY;
    use test_case::test_case;

    const LAB: &[usize] = &[4, 6, 2, 0, 8, 7, 5, 9, 3, 1];
    const PTN: &[usize] = &[
        NAUTY_INFINITY,
        3,
        NAUTY_INFINITY,
        1,
        2,
        NAUTY_INFINITY,
        NAUTY_INFINITY,
        0,
        2,
        0,
    ];
    #[test_case(LAB, PTN, 0, &[&[4, 6, 2, 0, 8, 7, 5, 9], &[3, 1]])]
    #[test_case(LAB, PTN, 1, &[&[4, 6, 2, 0], &[8, 7, 5, 9], &[3, 1]])]
    #[test_case(LAB, PTN, 2, &[&[4, 6, 2, 0], &[8], &[7, 5, 9], &[3], &[1]])]
    #[test_case(LAB, PTN, 3, &[&[4, 6], &[2, 0], &[8], &[7, 5, 9], &[3], &[1]])]
    fn test_partition(lab: &[usize], ptn: &[usize], level: usize, expected: &[&[usize]]) {
        assert_eq!(
            PartitionNest::new(Vec::from(lab), Vec::from(ptn)).partition(level),
            expected
        );
    }

    #[test_case(LAB, PTN, 3, "4, 6 | 2, 0 | 8 | 7, 5, 9 | 3 | 1")]
    fn test_partition_string(lab: &[usize], ptn: &[usize], level: usize, expected: &str) {
        assert_eq!(
            PartitionNest::new(Vec::from(lab), Vec::from(ptn)).partition_string(level),
            expected
        );
    }

    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", PartitionNest::new(Vec::from(LAB), Vec::from(PTN))),
            "".to_owned()
                + "PartitionNest { lab_ptn: [(4, 2000000002), (6, 3), (2, 2000000002), (0, 1), (8, 2), (7, 2000000002), (5, 2000000002), (9, 0), (3, 2), (1, 0)] }\n"
                + "0:  4, 6, 2, 0, 8, 7, 5, 9 | 3, 1\n"
                + "1:  4, 6, 2, 0 | 8, 7, 5, 9 | 3, 1\n"
                + "2:  4, 6, 2, 0 | 8 | 7, 5, 9 | 3 | 1\n"
                + "3:  4, 6 | 2, 0 | 8 | 7, 5, 9 | 3 | 1\n"
        );
    }
}
