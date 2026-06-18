use crate::nauty::NAUTY_INFINITY;
use std::{fmt::Debug, mem, ops::Index};

pub struct PartitionNest {
    lab: Vec<usize>,
    ptn: Vec<usize>,
    numcell: Vec<usize>,
    count: Vec<Vec<usize>>,
}

impl PartitionNest {
    pub fn new(lab: Vec<usize>, ptn: Vec<usize>) -> Self {
        assert_eq!(lab.len(), ptn.len());
        let mut nest = Self {
            lab,
            ptn,
            numcell: vec![],
            count: vec![],
        };
        let max_level = nest.max_level();
        if let Some(max_level) = max_level {
            for level in 0..=max_level {
                let partition = nest.partition_vec(level);
                nest.numcell.push(partition.len());
                nest.count.push(
                    nest.partition_vec(level)
                        .into_iter()
                        .map(|cell| cell.len())
                        .collect(),
                );
            }
        }
        nest
    }

    pub fn swap(&mut self, c1: usize, c2: usize) {
        self.lab.swap(c1, c2);
    }

    pub fn chunk_by(&self, level: usize) -> PartitionNestChunkBy {
        PartitionNestChunkBy {
            lab_slice: &self.lab,
            ptn_slice: &self.ptn,
            level,
        }
    }

    pub fn partition_vec(&self, level: usize) -> Vec<&[usize]> {
        self.chunk_by(level).collect::<Vec<_>>()
    }

    pub fn partition_string(&self, level: usize) -> String {
        itertools::Itertools::intersperse(
            self.partition_vec(level)
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

    pub fn max_level(&self) -> Option<usize> {
        self.ptn
            .iter()
            .filter(|l| l < &&NAUTY_INFINITY)
            .max()
            .copied()
    }

    /// last index of cell containing i for partition of given level
    pub fn cell_end(&self, i: usize, level: usize) -> usize {
        let mut end = i;
        while self.ptn[end] > level {
            end += 1;
        }
        end
    }

    pub fn get(&self, i: usize) -> usize {
        self.lab[i]
    }
}

impl Index<usize> for PartitionNest {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.lab[index]
    }
}

impl Debug for PartitionNest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PartitionNest")
            .field("lab", &self.lab)
            .field("ptn", &self.ptn)
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

pub struct Partition<'a> {
    nest: &'a PartitionNest,
    level: usize,
}

pub struct Cell<'a> {
    partition: &'a Partition<'a>,
    first_lab_index: usize,
}

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionNestChunkBy<'a> {
    lab_slice: &'a [usize],
    ptn_slice: &'a [usize],
    level: usize,
}

impl<'a> Iterator for PartitionNestChunkBy<'a> {
    type Item = &'a [usize];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.lab_slice.is_empty() {
            None
        } else {
            let mut len = 1;
            let ptn_iter = self.ptn_slice.iter();
            for ptn in ptn_iter {
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
        todo!()
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
            PartitionNest::new(Vec::from(lab), Vec::from(ptn)).partition_vec(level),
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

    #[test_case(LAB, PTN, 0, 0, 7)]
    #[test_case(LAB, PTN, 0, 1, 7)]
    #[test_case(LAB, PTN, 0, 7, 7)]
    #[test_case(LAB, PTN, 0, 8, 9)]
    #[test_case(LAB, PTN, 0, 9, 9)]
    #[test_case(LAB, PTN, 1, 0, 3)]
    #[test_case(LAB, PTN, 1, 3, 3)]
    #[test_case(LAB, PTN, 1, 4, 7)]
    #[test_case(LAB, PTN, 1, 7, 7)]
    #[test_case(LAB, PTN, 1, 8, 9)]
    #[test_case(LAB, PTN, 1, 9, 9)]
    #[test_case(LAB, PTN, 2, 0, 3)]
    #[test_case(LAB, PTN, 2, 3, 3)]
    #[test_case(LAB, PTN, 2, 4, 4)]
    #[test_case(LAB, PTN, 2, 5, 7)]
    #[test_case(LAB, PTN, 2, 7, 7)]
    #[test_case(LAB, PTN, 2, 8, 8)]
    #[test_case(LAB, PTN, 2, 9, 9)]
    #[test_case(LAB, PTN, 3, 0, 1)]
    #[test_case(LAB, PTN, 3, 1, 1)]
    #[test_case(LAB, PTN, 3, 2, 3)]
    #[test_case(LAB, PTN, 3, 3, 3)]
    #[test_case(LAB, PTN, 3, 4, 4)]
    #[test_case(LAB, PTN, 3, 5, 7)]
    #[test_case(LAB, PTN, 3, 7, 7)]
    #[test_case(LAB, PTN, 3, 8, 8)]
    #[test_case(LAB, PTN, 3, 9, 9)]
    fn test_split2(lab: &[usize], ptn: &[usize], level: usize, split1: usize, split2: usize) {
        assert_eq!(
            PartitionNest::new(Vec::from(lab), Vec::from(ptn)).cell_end(split1, level),
            split2
        );
    }

    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", PartitionNest::new(Vec::from(LAB), Vec::from(PTN))),
            "".to_owned()
                + "PartitionNest { lab: [4, 6, 2, 0, 8, 7, 5, 9, 3, 1], ptn: [2000000002, 3, 2000000002, 1, 2, 2000000002, 2000000002, 0, 2, 0] }\n"
                + "0:  4, 6, 2, 0, 8, 7, 5, 9 | 3, 1\n"
                + "1:  4, 6, 2, 0 | 8, 7, 5, 9 | 3, 1\n"
                + "2:  4, 6, 2, 0 | 8 | 7, 5, 9 | 3 | 1\n"
                + "3:  4, 6 | 2, 0 | 8 | 7, 5, 9 | 3 | 1\n"
        );
    }
}
