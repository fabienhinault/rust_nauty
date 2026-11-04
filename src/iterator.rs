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

/// borrowed from partition_in_place for rust,
/// https://doc.rust-lang.org/src/core/iter/traits/iterator.rs.html#2334-2337
/// but tries to mimic the algorithm used in nauty's code
pub(crate) fn partition_in_place_orig<P>(slice: &mut [usize], mut predicate: P) -> usize
where
    P: FnMut(&usize) -> bool,
{
    let mut j0 = 0;
    let mut j1 = slice.len() - 1;
    while j0 <= j1 {
        if predicate(&slice[j0]) {
            j0 += 1;
        } else {
            slice.swap(j0, j1);
            if j1 == 0 {
                break;
            }
            j1 -= 1;
        }
    }
    j0
}
