use crate::{
    nautil::SetWordNautilTrait,
    nauty::{Graph, SetTrait},
};

/*****************************************************************************
*                                                                            *
*  isautom(g,perm,digraph,m,n) = TRUE iff perm is an automorphism of g       *
*  (i.e., g^perm = g).  Symmetry is assumed unless digraph = TRUE.           *
*                                                                            *
*****************************************************************************/
fn is_autom(g: &Graph, perm: &[usize]) -> bool {
    for (i_row, row) in g.0.iter().enumerate() {
        let p_row = g.get(perm[i_row]);
        for pos in row.masked(i_row).ones_iter() {
            if !p_row[perm[pos]] {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::nauty::test::create_n_circle;

    #[test]
    fn test_is_autom() {
        assert!(is_autom(&create_n_circle(4), &[1, 2, 3, 0]));
        assert!(!is_autom(&create_n_circle(4), &[1, 0, 2, 3]));
    }
}
