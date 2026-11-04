use crate::{
    iterator::partition_in_place,
    nauty::{
        Set, VecMap,
        graph::{self, Graph},
        partition_nest::partition::cell::Cell,
    },
};
use std::{
    iter::once,
    mem::swap,
    ops::{Index, IndexMut},
};
/// struct of results of a given function on a usize range.
/// Used internally for cell split.
struct FResults {
    f_results: Vec<usize>,
    f_antecedant_counts: VecMap,
    min: usize,
    max: usize,
}

impl FResults {
    fn is_const(&self) -> Option<usize> {
        (self.min == self.max).then_some(self.min)
    }
}

pub struct CellMut<'a> {
    pub(crate) level: usize,
    pub(crate) first_lab_index: usize,
    pub(crate) cell_lab: &'a mut [usize],
    pub(crate) cell_ptn: &'a mut [usize],
}

impl<'a> CellMut<'a> {
    pub fn len(&self) -> usize {
        self.cell_lab.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.cell_lab.is_empty()
    }

    pub fn to_vec(&self) -> Vec<usize> {
        self.cell_lab.to_vec()
    }

    pub fn swap(&mut self, a: usize, b: usize) {
        self.cell_lab.swap(a, b)
    }

    pub fn raw_split(&mut self, i: usize) {
        self.cell_ptn[i] = self.level;
    }

    pub fn partition_index(&self, i: usize) -> usize {
        self.first_lab_index + i
    }

    pub fn copy_from_slice(&mut self, src: &[usize]) {
        self.cell_lab.copy_from_slice(src)
    }

    /// partition the cell, put splitter neighbours first
    /// returns the number of neighbours
    pub fn partition_in_place_std(&mut self, splitter_neighbours: &Set) -> usize {
        partition_in_place(self.cell_lab, |c| splitter_neighbours[*c])
    }

    /// partition the cell, put splitter neighbours first
    /// returns the number of neighbours
    pub fn partition_in_place_orig(&mut self, splitter_neighbours: &Set) -> usize {
        let mut c1: usize = 0;
        // c2 can be -1
        let mut c2: isize = self.len() as isize - 1;
        while c1 as isize <= c2 {
            if splitter_neighbours[self[c1]] {
                c1 += 1;
            } else {
                self.swap(c1, c2 as usize);
                c2 -= 1;
            }
        }
        c1
    }

    pub fn split_trivial(&mut self, splitter_neighbours: &Set) -> usize {
        let i = self.partition_in_place_orig(splitter_neighbours);
        if i > 0 && i < self.len() {
            self.raw_split(i - 1);
        }
        i
    }

    pub fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.cell_lab.iter()
    }

    pub fn split_from_fn<F: FnMut(&usize) -> usize>(&mut self, f: F) -> SplitResult {
        let f_results = self.get_f_results(f);
        if let Some(value) = f_results.is_const() {
            return SplitResult::Const(value);
        }
        self.split_from_f_results(f_results)
    }

    fn get_f_results<F: FnMut(&usize) -> usize>(&self, f: F) -> FResults {
        let f_results: Vec<usize> = self.iter().map(f).collect();
        let (min, max, f_antecedant_counts) = f_results
            .iter()
            .fold((None, None, VecMap::new()), fold_step);
        let min = min.expect("min");
        let max = max.expect("max");
        FResults {
            f_results,
            f_antecedant_counts,
            min,
            max,
        }
    }

    fn split_from_f_results(&mut self, f_results: FResults) -> SplitResult {
        let FResults {
            f_results,
            f_antecedant_counts,
            min,
            max,
        } = f_results;
        let mut f_value_indices = VecMap::new();
        let mut current_f_value_index: usize = 0;
        let mut biggest_cell_size = 0;
        let mut biggest_cell_pos = None;
        for f_value in min..=max {
            if let Some(count) = f_antecedant_counts.get_safely(f_value) {
                if count > biggest_cell_size {
                    biggest_cell_size = count;
                    biggest_cell_pos = Some(current_f_value_index);
                }
                f_value_indices.set(f_value, current_f_value_index);
                current_f_value_index += count;
                if current_f_value_index < self.len() {
                    self.cell_ptn[current_f_value_index - 1] = self.level;
                }
            }
        }
        let result = SplitResult::Split {
            biggest_cell_pos: biggest_cell_pos.expect("biggest_cell_pos"),
            cells_indices: f_value_indices
                .0
                .iter()
                .flatten()
                .copied()
                .chain(once(self.len()))
                .collect(),
            cells_values: f_value_indices
                .0
                .iter()
                .enumerate()
                .filter_map(|(i, o)| o.map(|_| i))
                .collect(),
        };
        let mut new_lab = vec![0; self.len()];
        for (result, lab) in f_results.iter().zip(self.iter()) {
            new_lab[f_value_indices.get(*result)] = *lab;
            f_value_indices.increment(*result);
        }
        self.copy_from_slice(&new_lab);
        result
    }

    fn splitters_count(&self, splitters: &Set, g: &graph::Graph, i: usize) -> usize {
        (splitters.clone() & g.0[i].clone()).count_ones()
    }

    pub fn split_from_splitters(&mut self, splitters: &Set, g: &graph::Graph) -> SplitResult {
        let f_results = self.get_f_results(|i: &usize| self.splitters_count(splitters, g, *i));
        if let Some(value) = f_results.is_const() {
            return SplitResult::Const(value);
        }
        self.split_from_f_results(f_results)
    }

    /*****************************************************************************
     *                                                                            *
     *  breakout(lab,ptn,level,tc,tv,active,m) operates on the partition at       *
     *  the specified level in the partition nest (lab,ptn).  It finds the        *
     *  element tv, which is in the cell C starting at index tc in lab (it had    *
     *  better be) and splits C in the two cells {tv} and C\{tv}, in that order.  *
     *  It also sets the set active to contain just the element tc.               *
     *                                                                            *
     *  GLOBALS ACCESSED: bit<r>                                                  *
     *                                                                            *
     *****************************************************************************/
    // nautil.c 613
    /// Operate at level + 1, isolate vertex tv in one cell, and other ones
    /// in another.
    pub fn breakout(&mut self, tv: usize) {
        let mut tmp = tv;
        swap(&mut tmp, &mut self[0]);
        let mut i = 1;
        while tmp != tv {
            swap(&mut tmp, &mut self[i]);
            i += 1;
        }
        self.cell_ptn[0] = self.level + 1;
    }
}

pub enum SplitResult {
    Const(usize),
    Split {
        biggest_cell_pos: usize,
        cells_indices: Vec<usize>,
        cells_values: Vec<usize>,
    },
}

fn fold_step(
    acc: (Option<usize>, Option<usize>, VecMap),
    &x: &usize,
) -> (Option<usize>, Option<usize>, VecMap) {
    let (min, max, mut f_antecedant_counts) = acc;
    let c = f_antecedant_counts.get_safely(x).unwrap_or(0);
    f_antecedant_counts.set(x, c + 1);
    (
        min.map_or(Some(x), |m: usize| Some(m.min(x))),
        max.map_or(Some(x), |m: usize| Some(m.max(x))),
        f_antecedant_counts,
    )
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

#[cfg(test)]
mod test {
    use super::*;
    use crate::nauty::NAUTY_INFINITY;
    use crate::nauty::partition_nest::PartitionNest;
    use crate::nauty::test::create_from_offsets;

    #[test]
    fn test_split_from_fn() {
        let nest = PartitionNest::new(
            vec![0, 3, 4, 2, 1, 6, 5],
            vec![
                2,
                NAUTY_INFINITY,
                2,
                NAUTY_INFINITY,
                NAUTY_INFINITY,
                NAUTY_INFINITY,
                0,
            ],
        );
        let mut partition = nest.partition(2);
        assert_eq!(partition.numcells(), 3);
        let mut cells_mut = partition.cells_mut();
        let mut cell_mut = cells_mut.nth(2).unwrap();
        let f = |vertex: &usize| if [1, 6].contains(vertex) { 1 } else { 0 };
        cell_mut.split_from_fn(f);
        assert_eq!(partition.numcells(), 4);
        assert_eq!(
            partition.nest.clone(),
            PartitionNest::new(
                vec![0, 3, 4, 2, 5, 1, 6],
                vec![2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 0],
            )
        );
    }

    #[test]
    fn test_split_from_cell() {
        let nest = PartitionNest::new(
            vec![0, 3, 4, 2, 1, 6, 5],
            vec![
                2,
                NAUTY_INFINITY,
                2,
                NAUTY_INFINITY,
                NAUTY_INFINITY,
                NAUTY_INFINITY,
                0,
            ],
        );
        let mut partition = nest.partition(2);
        assert_eq!(partition.numcells(), 3);
        let mut cells = partition.cells();
        let cell = cells.nth(1).unwrap();
        let splitters = cell.get_splitters();
        let mut cells_mut = partition.cells_mut();
        let mut cell_mut = cells_mut.nth(2).unwrap();
        let g = create_from_offsets(7, &[3, 4]);
        cell_mut.split_from_splitters(&splitters, &g);
        assert_eq!(partition.numcells(), 4);
        assert_eq!(
            partition.nest.clone(),
            PartitionNest::new(
                vec![0, 3, 4, 2, 5, 1, 6],
                vec![2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 2, NAUTY_INFINITY, 0],
            )
        );
    }

    #[test]
    fn test_breakout() {
        let nest = PartitionNest::new(vec![2, 1, 0], vec![NAUTY_INFINITY, NAUTY_INFINITY, 0]);
        let mut partition = nest.partition(1);
        partition.get_cell_mut(0).breakout(0);
        assert_eq!(
            partition.nest,
            PartitionNest::new(vec![0, 2, 1], vec![2, NAUTY_INFINITY, 0])
        );
    }
}
