use crate::nauty::{NAUTY_INFINITY, Set};
use partition_nest_chunk_by::PartitionNestChunkBy;
use partition_nest_chunk_by_mut::PartitionNestChunkByMut;
use std::{fmt::Debug, ops::Index};

pub mod partition;
mod partition_nest_chunk_by;
mod partition_nest_chunk_by_mut;

/// *    A partition nest is represented by a pair (lab,ptn), where lab and ptn  *
/// *    are int arrays.  The "partition at level x" is the partition whose      *
/// *    cells are {lab[i],lab[i+1],...,lab[j]}, where [i,j] is a maximal        *
/// *    subinterval of [0,n-1] such that ptn[k] > x for i <= k < j and          *
/// *    ptn[j] <= x.  The partition at level 0 is given to nauty by the user.   *
/// *    This is  refined for the root of the tree, which has level 1.           *
#[derive(Default, PartialEq)]
pub struct PartitionNest {
    /// lab must always be a permutation of [[0, n-1]]
    lab: Vec<usize>,
    ptn: Vec<usize>,
    /// Denormalized numbers of cells at all levels.
    numcells: Vec<usize>,
}

impl PartitionNest {
    pub fn new(lab: Vec<usize>, ptn: Vec<usize>) -> Self {
        assert_eq!(lab.len(), ptn.len());
        let mut nest = Self {
            lab,
            ptn,
            numcells: vec![],
        };
        for level in 0..=nest.max_level() {
            let partition = nest.partition_vec(level);
            nest.numcells.push(partition.len());
        }
        nest
    }

    pub fn assert_is_sane(&self) {
        let mut lab = self.lab.clone();
        lab.sort();
        assert_eq!(lab, (0..lab.len()).collect::<Vec<_>>());
    }

    pub fn swap(&mut self, c1: usize, c2: usize) {
        self.lab.swap(c1, c2);
    }

    pub fn chunk_by(&self, level: usize) -> PartitionNestChunkBy<'_> {
        PartitionNestChunkBy::new(&self.lab, &self.ptn, level)
    }

    pub fn chunk_by_mut(&mut self, level: usize) -> PartitionNestChunkByMut<'_> {
        PartitionNestChunkByMut::new(&mut self.lab, &mut self.ptn, level)
    }

    pub fn partition(self, level: usize) -> partition::Partition {
        partition::Partition { nest: self, level }
    }

    pub fn partition_vec(&self, level: usize) -> Vec<&[usize]> {
        self.chunk_by(level).collect::<Vec<_>>()
    }

    pub fn partition_string(&self, level: usize) -> String {
        itertools::Itertools::intersperse(
            self.chunk_by(level)
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

    pub fn max_level(&self) -> usize {
        *self
            .ptn
            .iter()
            .filter(|l| l < &&NAUTY_INFINITY)
            .max()
            .unwrap_or(&0)
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

    pub fn numcells(&self, level: usize) -> usize {
        self.numcells[level.min(self.max_level())]
    }

    pub fn max_level_numcells(&self) -> usize {
        self.numcells[self.max_level()]
    }

    fn numcells_mut(&mut self, level: usize) -> &mut usize {
        self.extend_numcells(level);
        &mut self.numcells[level]
    }

    fn extend_numcells(&mut self, level: usize) {
        if level >= self.numcells.len() {
            self.numcells.extend_from_slice(&vec![
                self.max_level_numcells();
                level - self.numcells.len() + 1
            ]);
        }
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
        writeln!(f)?;
        for l in 0..=self.max_level() {
            writeln!(f, "{l}:  {}", self.partition_string(l))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::nauty::{NAUTY_INFINITY, partition_nest::partition::Partition};
    use bitvec::{bitvec, order::Msb0};
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

    #[test]
    fn test_unpartioned_4() {
        assert_eq!(
            format!(
                "{:?}",
                PartitionNest::new(
                    vec![0, 1, 2, 3],
                    vec![NAUTY_INFINITY, NAUTY_INFINITY, NAUTY_INFINITY, 0]
                )
            ),
            "".to_owned()
                + "PartitionNest { lab: [0, 1, 2, 3], ptn: [2000000002, 2000000002, 2000000002, 0] }\n"
                + "0:  0, 1, 2, 3\n"
        );
    }

    #[test_case(Partition::new(PartitionNest::new(vec![0, 2, 1], vec![2, NAUTY_INFINITY, 0]),2), 2, 1, bitvec![usize, Msb0; 0; 3], 0, &[0, 1, 2], &[2, NAUTY_INFINITY, 0],2)]
    #[allow(clippy::too_many_arguments)]
    fn test_partition_inplace_orig(
        mut partition: partition::Partition,
        expected_numcells_before: usize,
        i_cell: usize,
        gptr: Set,
        expected_true_count: usize,
        expected_lab: &[usize],
        expected_ptn: &[usize],
        expected_numcells_after: usize,
    ) {
        assert_eq!(partition.numcells(), expected_numcells_before);
        let mut cells_mut = partition.cells_mut();
        let mut cell = cells_mut.next().unwrap();
        for _ in 0..(i_cell) {
            cell = cells_mut.next().unwrap();
        }
        let true_count = cell.partition_in_place_orig(&gptr);
        assert_eq!(true_count, expected_true_count);
        assert_eq!(partition.numcells(), expected_numcells_after);
        assert_eq!(partition.nest.lab, expected_lab);
        assert_eq!(partition.nest.ptn, expected_ptn);
    }
}
