use crate::nauty::{Set, VecMap};
use std::ops::{Index, IndexMut};

/// https://doc.rust-lang.org/src/core/iter/traits/iterator.rs.html#2334-2337
///
/// Reorders the elements of this iterator *in-place* according to the given predicate,
/// such that all those that return `true` precede all those that return `false`.
/// Returns the number of `true` elements found.
///
/// The relative order of partitioned items is not maintained.
///
/// # Current implementation
///
/// The current algorithm tries to find the first element for which the predicate evaluates
/// to false and the last element for which it evaluates to true, and repeatedly swaps them.
///
/// Time complexity: *O*(*n*)
///
/// See also [`is_partitioned()`] and [`partition()`].
///
/// [`is_partitioned()`]: Iterator::is_partitioned
/// [`partition()`]: Iterator::partition
///
/// # Examples
///
/// ```
/// #![feature(iter_partition_in_place)]
///
/// let mut a = [1, 2, 3, 4, 5, 6, 7];
///
/// // Partition in-place between evens and odds
/// let i = a.iter_mut().partition_in_place(|n| n % 2 == 0);
///
/// assert_eq!(i, 3);
/// assert!(a[..i].iter().all(|n| n % 2 == 0)); // evens
/// assert!(a[i..].iter().all(|n| n % 2 == 1)); // odds
/// ```
pub(crate) fn partition_in_place<P>(slice: &mut [usize], ref mut predicate: P) -> usize
where
    P: FnMut(&usize) -> bool,
{
    // FIXME: should we worry about the count overflowing? The only way to have more than
    // `usize::MAX` mutable references is with ZSTs, which aren't useful to partition...

    // These closure "factory" functions exist to avoid genericity in `Self`.

    #[inline]
    fn is_false<'a, T>(
        predicate: &'a mut impl FnMut(&T) -> bool,
        true_count: &'a mut usize,
    ) -> impl FnMut(&&mut T) -> bool + 'a {
        move |x| {
            let p = predicate(&**x);
            *true_count += p as usize;
            !p
        }
    }

    #[inline]
    fn is_true<T>(predicate: &mut impl FnMut(&T) -> bool) -> impl FnMut(&&mut T) -> bool + '_ {
        move |x| predicate(&**x)
    }

    let mut iter = slice.iter_mut();
    // Repeatedly find the first `false` and swap it with the last `true`.
    let mut true_count = 0;
    while let Some(head) = iter.find(is_false(predicate, &mut true_count)) {
        if let Some(tail) = iter.rfind(is_true(predicate)) {
            std::mem::swap(head, tail);
            true_count += 1;
        } else {
            break;
        }
    }
    true_count
}

pub struct CellMut<'a> {
    pub(crate) level: usize,
    pub(crate) first_lab_index: usize,
    pub(crate) cell_lab: &'a mut [usize],
    pub(crate) cell_ptn: &'a mut [usize],
    pub(crate) numcells: &'a mut usize,
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

    pub fn split_trivial(&mut self, i: usize) {
        self.cell_ptn[i] = self.level;
        *self.numcells += 1;
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

    pub fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.cell_lab.iter()
    }

    pub fn split_from_fn<F: FnMut(&usize) -> usize>(&mut self, f: F) {
        let (min, max, f_results, f_antecedant_counts) = self
            .iter()
            .map(f)
            .fold((None, None, vec![], VecMap::new()), fold_step);
        let min = min.expect("min");
        let max = max.expect("max");
        if min == max {
            return;
        }
        let mut f_value_indices = VecMap::new();
        let mut current_f_value_index: usize = 0;
        for f_value in min..=max {
            if let Some(count) = f_antecedant_counts.get_safely(f_value) {
                f_value_indices.set(f_value, current_f_value_index);
                if current_f_value_index != 0 {
                    *self.numcells += 1;
                }
                current_f_value_index += count;
                if current_f_value_index < self.len() {
                    self.cell_ptn[current_f_value_index - 1] = self.level;
                }
            }
        }
        let mut new_lab = vec![0; self.len()];
        for (result, lab) in f_results.iter().zip(self.iter()) {
            new_lab[f_value_indices.get(*result)] = *lab;
            f_value_indices.increment(*result);
        }
        self.copy_from_slice(&new_lab);
    }
}

fn fold_step(
    acc: (Option<usize>, Option<usize>, Vec<usize>, VecMap),
    x: usize,
) -> (Option<usize>, Option<usize>, Vec<usize>, VecMap) {
    let (min, max, mut f_results, mut f_antecedant_counts) = acc;
    f_results.push(x);
    let c = f_antecedant_counts.get_safely(x).unwrap_or(0);
    f_antecedant_counts.set(x, c + 1);
    (
        min.map_or(Some(x), |m: usize| Some(m.min(x))),
        max.map_or(Some(x), |m: usize| Some(m.max(x))),
        f_results,
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

    #[test]
    fn test() {
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
}
