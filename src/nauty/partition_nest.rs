use crate::nauty::NAUTY_INFINITY;
use partition_nest_chunk_by::PartitionNestChunkBy;
use partition_nest_chunk_by_mut::PartitionNestChunkByMut;
use std::{
    fmt::Debug,
    ops::{Index, IndexMut},
};

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

    pub fn partition(self, level: usize) -> Partition {
        Partition { nest: self, level }
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

#[derive(Default, PartialEq, Debug)]
pub struct Partition {
    nest: PartitionNest,
    level: usize,
}

impl Partition {
    pub fn new(mut nest: PartitionNest, level: usize) -> Self {
        nest.extend_numcells(level);
        Self { nest, level }
    }

    pub fn numcells(&self) -> usize {
        self.nest.numcells(self.level)
    }

    pub fn numcells_mut(&mut self) -> &mut usize {
        self.nest.numcells_mut(self.level)
    }

    /// last index of cell containing i for partition of given level
    pub fn cell_end(&self, i: usize) -> usize {
        let mut end = i;
        while self.nest.ptn[end] > self.level {
            end += 1;
        }
        end
    }

    pub fn raw_cells(&self) -> PartitionNestChunkBy<'_> {
        self.nest.chunk_by(self.level)
    }

    pub fn cells(&self) -> PartitionCellsIter<'_> {
        PartitionCellsIter {
            partition: self,
            index: 0,
        }
    }

    pub fn cells_mut(&mut self) -> PartitionCellsIterMut<'_> {
        PartitionCellsIterMut {
            partition: self,
            index: 0,
        }
    }

    pub fn split(&mut self, i: usize) {
        self.nest.lab[i] = self.level;
        self.nest.numcells[self.level] += 1;
    }

    fn len(&self) -> usize {
        self.nest.lab.len()
    }
}

impl Index<usize> for Partition {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.nest.lab[index]
    }
}

impl Default for &mut Partition {
    fn default() -> Self {
        todo!()
    }
}

pub struct Cell<'a> {
    partition: &'a Partition,
    first_lab_index: usize,
    cell_lab: &'a [usize],
}

impl<'a> Cell<'a> {
    pub fn is_at_end(&self) -> bool {
        self.first_lab_index == self.partition.nest.lab.len()
    }

    pub fn len(&self) -> usize {
        self.cell_lab.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.cell_lab.is_empty()
    }
}

pub struct CellMut<'a> {
    level: usize,
    first_lab_index: usize,
    cell_lab: &'a mut [usize],
    cell_ptn: &'a mut [usize],
    numcells: &'a mut usize,
}

impl<'a> CellMut<'a> {
    pub fn len(&self) -> usize {
        self.cell_lab.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.cell_lab.is_empty()
    }

    pub fn swap(&mut self, a: usize, b: usize) {
        self.cell_lab.swap(a, b)
    }

    pub fn split(&mut self, i: usize) {
        self.cell_ptn[i] = self.level;
        *self.numcells += 1;
    }

    pub fn partition_index(&self, i: usize) -> usize {
        self.first_lab_index + i
    }
}

impl<'a> Index<usize> for CellMut<'a> {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.cell_lab[index]
    }
}

impl<'a> IndexMut<usize> for CellMut<'a> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.cell_lab[index]
    }
}

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionCellsIter<'a> {
    partition: &'a Partition,
    index: usize,
}

impl<'a> Iterator for PartitionCellsIter<'a> {
    type Item = Cell<'a>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.index == self.partition.nest.lab.len() {
            None
        } else {
            let mut len = 1;
            let ptn_iter = self.partition.nest.ptn[self.index..].iter();
            for ptn in ptn_iter {
                if *ptn > self.partition.level {
                    len += 1
                } else {
                    break;
                }
            }
            let first_lab_index = self.index;
            self.index += len;
            Some(Cell {
                partition: self.partition,
                first_lab_index,
                cell_lab: &self.partition.nest.lab[first_lab_index..self.index],
            })
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.index == self.partition.nest.lab.len() {
            (0, Some(0))
        } else {
            (1, Some(self.partition.nest.lab.len() - self.index))
        }
    }

    #[inline]
    fn last(self) -> Option<Self::Item> {
        todo!()
    }
}

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionCellsIterMut<'a> {
    partition: &'a mut Partition,
    index: usize,
}

// do not implement trait Iterator because of lifetime issue
impl<'a> PartitionCellsIterMut<'a> {
    pub fn next<'b>(&'b mut self) -> Option<CellMut<'b>> {
        if self.index == self.partition.len() {
            None
        } else {
            let mut len = 1;
            for ptn in self.partition.nest.ptn[self.index..].iter() {
                if *ptn > self.partition.level {
                    len += 1
                } else {
                    break;
                }
            }
            let first_lab_index = self.index;
            self.index += len;
            let lab_slice: &'b mut [usize] =
                &mut self.partition.nest.lab[first_lab_index..self.index];
            let ptn_slice: &'b mut [usize] =
                &mut self.partition.nest.ptn[first_lab_index..self.index];
            Some(CellMut {
                level: self.partition.level,
                first_lab_index,
                cell_lab: lab_slice,
                cell_ptn: ptn_slice,
                numcells: &mut self.partition.nest.numcells[self.partition.level],
            })
        }
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
}
