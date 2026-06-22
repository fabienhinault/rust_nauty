use crate::nauty::NAUTY_INFINITY;
use partition_nest_chunk_by::PartitionNestChunkBy;
use partition_nest_chunk_by_mut::PartitionNestChunkByMut;
use std::{fmt::Debug, mem, ops::Index};

mod partition_nest_chunk_by;
mod partition_nest_chunk_by_mut;

/// *    A partition nest is represented by a pair (lab,ptn), where lab and ptn  *
/// *    are int arrays.  The "partition at level x" is the partition whose      *
/// *    cells are {lab[i],lab[i+1],...,lab[j]}, where [i,j] is a maximal        *
/// *    subinterval of [0,n-1] such that ptn[k] > x for i <= k < j and          *
/// *    ptn[j] <= x.  The partition at level 0 is given to nauty by the user.   *
/// *    This is  refined for the root of the tree, which has level 1.           *
#[derive(Default)]
pub struct PartitionNest {
    lab: Vec<usize>,
    ptn: Vec<usize>,
    /// Denormalized numbers of cells at all levels.
    numcells: Vec<usize>,
    /// Denormalized counts of vertices in cells at all levels.
    /// count[level][i_cell] is the vertex count of cell i_cell at given level.
    count: Vec<Vec<usize>>,
}

impl PartitionNest {
    pub fn new(lab: Vec<usize>, ptn: Vec<usize>) -> Self {
        assert_eq!(lab.len(), ptn.len());
        let mut nest = Self {
            lab,
            ptn,
            numcells: vec![],
            count: vec![],
        };
        let max_level = nest.max_level();
        if let Some(max_level) = max_level {
            for level in 0..=max_level {
                let partition = nest.partition_vec(level);
                nest.numcells.push(partition.len());
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
        PartitionNestChunkBy::new(&self.lab, &self.ptn, level)
    }

    pub fn chunk_by_mut(&mut self, level: usize) -> PartitionNestChunkByMut {
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

#[derive(Default)]
pub struct Partition {
    nest: PartitionNest,
    level: usize,
}

impl Partition {
    pub fn numcells(&self) -> usize {
        self.nest.numcells[self.level]
    }

    /// last index of cell containing i for partition of given level
    pub fn cell_end(&self, i: usize) -> usize {
        let mut end = i;
        while self.nest.ptn[end] > self.level {
            end += 1;
        }
        end
    }

    pub fn raw_cells(&self) -> PartitionNestChunkBy {
        self.nest.chunk_by(self.level)
    }

    pub fn cells(&self) -> PartitionCellsIter {
        PartitionCellsIter {
            partition: self,
            index: 0,
        }
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
    cell_ptn: &'a [usize],
}

impl<'a> Cell<'a> {
    pub fn is_at_end(&self) -> bool {
        self.first_lab_index == self.partition.nest.lab.len()
    }

    pub fn len(&self) -> usize {
        self.cell_lab.len()
    }
}

pub struct CellMut<'a> {
    level:
    first_lab_index: usize,
    cell_lab: &'a mut [usize],
    cell_ptn: &'a mut [usize],
}

pub struct CellMove {
    partition: Partition,
    first_lab_index: usize,
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
            self.index = self.index + len;
            Some(Cell {
                partition: self.partition.clone(),
                first_lab_index,
                cell_lab: &self.partition.nest.lab[first_lab_index..self.index],
                cell_ptn: &self.partition.nest.ptn[first_lab_index..self.index],
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
    fn last(mut self) -> Option<Self::Item> {
        todo!()
    }
}

#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct PartitionCellsIterMut<'a> {
    partition: &'a mut Partition,
    index: usize,
}

impl<'a> Iterator for PartitionCellsIterMut<'a> {
    type Item = CellMut<'a>;

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
            self.index = self.index + len;
            let partition: &mut Partition = mem::take(&mut self.partition);
            Some(CellMut {
                first_lab_index,
                cell_lab: &mut partition.nest.lab[first_lab_index..self.index],
                cell_ptn: &mut partition.nest.ptn[first_lab_index..self.index],
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
