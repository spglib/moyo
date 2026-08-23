use std::collections::BTreeMap;

use serde::Serialize;

use crate::base::{MoyoError, Operations};
use crate::subgroup::finite_group::FiniteGroup;

/// A Translationengleiche subgroup embedded in its parent space group.
///
/// A Translationengleiche subgroup has the same translation lattice as its
/// parent. Consequently, it is the inverse image of a subgroup of the parent
/// point group under the projection from the space group to its point group.
#[derive(Debug, Clone, Serialize)]
pub struct TranslationengleicheSubgroup {
    /// Coset representatives of the subgroup, in the parent primitive basis.
    pub operations: Operations,
    /// Mapping from [`operations`](Self::operations) to the input parent operations.
    pub operation_indices: Vec<usize>,
    /// Translationengleiche index `[G/T : H/T]`.
    pub translationengleiche_index: usize,
    /// Shortest covering-chain length from the parent in the subgroup lattice.
    ///
    /// The parent itself has depth zero and its maximal proper
    /// Translationengleiche subgroups have depth one.
    pub depth: usize,
}

/// A member of a parent-conjugacy class of embedded Translationengleiche subgroups.
#[derive(Debug, Clone, Serialize)]
pub struct TranslationengleicheSubgroupConjugate {
    /// The conjugate subgroup embedded in the input parent.
    pub subgroup: TranslationengleicheSubgroup,
    /// Index of an input parent operation whose inverse conjugates the
    /// representative to [`subgroup`](Self::subgroup).
    pub conjugator_index: usize,
}

/// A parent-conjugacy class of embedded Translationengleiche subgroups.
#[derive(Debug, Clone, Serialize)]
pub struct TranslationengleicheSubgroupConjugacyClass {
    /// Deterministically selected representative of this conjugacy class.
    pub representative: TranslationengleicheSubgroup,
    /// Every distinct conjugate embedding, including the representative first.
    pub conjugates: Vec<TranslationengleicheSubgroupConjugate>,
    /// Indices of input parent operations that normalize the representative.
    pub normalizer_indices: Vec<usize>,
    /// Permutations induced on the representative operations by the normalizer.
    ///
    /// This is aligned with [`normalizer_indices`](Self::normalizer_indices).
    /// For a permutation `p`, operation position `i` is mapped to `p[i]` by
    /// inverse conjugation with the corresponding normalizer operation.
    pub normalizer_permutations: Vec<Vec<usize>>,
}

/// Enumerate Translationengleiche subgroups of a primitive space group.
///
/// `prim_operations` must contain exactly one representative for every point-
/// group element. The operations are validated to form a group modulo integer
/// lattice translations using `epsilon`.
///
/// The result contains all Translationengleiche subgroups, including the
/// parent and its pure-translation subgroup, partitioned into conjugacy classes
/// under the parent space group. Each subgroup is returned in the input
/// primitive basis and its operation indices refer directly to
/// `prim_operations`.
pub fn enumerate_translationengleiche_subgroups(
    prim_operations: &Operations,
    epsilon: f64,
) -> Result<Vec<TranslationengleicheSubgroupConjugacyClass>, MoyoError> {
    let finite_group = FiniteGroup::from_operations(prim_operations, epsilon)?;
    let subgroups = finite_group.enumerate_subgroups();
    let depths = subgroup_depths(&subgroups);

    let mut classes = Vec::new();
    for conjugacy_class in finite_group.conjugacy_classes(&subgroups) {
        let representative_indices = &conjugacy_class.representative;
        let representative = make_translationengleiche_subgroup(
            representative_indices,
            prim_operations,
            finite_group.order(),
            &depths,
        );
        let conjugates = conjugacy_class
            .conjugates
            .into_iter()
            .map(|conjugate| TranslationengleicheSubgroupConjugate {
                subgroup: make_translationengleiche_subgroup(
                    &conjugate.subgroup,
                    prim_operations,
                    finite_group.order(),
                    &depths,
                ),
                conjugator_index: conjugate.conjugator_index,
            })
            .collect::<Vec<_>>();
        classes.push(TranslationengleicheSubgroupConjugacyClass {
            representative,
            conjugates,
            normalizer_indices: conjugacy_class.normalizer_indices,
            normalizer_permutations: conjugacy_class.normalizer_permutations,
        });
    }

    Ok(classes)
}

fn make_translationengleiche_subgroup(
    indices: &[usize],
    prim_operations: &Operations,
    parent_order: usize,
    depths: &BTreeMap<Vec<usize>, usize>,
) -> TranslationengleicheSubgroup {
    TranslationengleicheSubgroup {
        operations: indices
            .iter()
            .map(|&index| prim_operations[index].clone())
            .collect(),
        operation_indices: indices.to_vec(),
        translationengleiche_index: parent_order / indices.len(),
        depth: depths[indices],
    }
}

fn subgroup_depths(subgroups: &[Vec<usize>]) -> BTreeMap<Vec<usize>, usize> {
    let mut ordered_subgroups = subgroups.iter().collect::<Vec<_>>();
    ordered_subgroups.sort_by(|lhs, rhs| rhs.len().cmp(&lhs.len()).then_with(|| lhs.cmp(rhs)));
    let Some((parent, proper_subgroups)) = ordered_subgroups.split_first() else {
        return BTreeMap::new();
    };

    let mut depths = BTreeMap::new();
    depths.insert((*parent).clone(), 0);

    for &subgroup in proper_subgroups {
        let covering_supergroups = ordered_subgroups.iter().copied().filter(|candidate| {
            is_proper_subset(subgroup, candidate)
                && !ordered_subgroups.iter().copied().any(|between| {
                    is_proper_subset(subgroup, between) && is_proper_subset(between, candidate)
                })
        });
        let depth = covering_supergroups
            .map(|parent| depths[parent] + 1)
            .min()
            .expect("every proper subgroup has a covering supergroup");
        depths.insert(subgroup.clone(), depth);
    }

    depths
}

fn is_proper_subset(lhs: &[usize], rhs: &[usize]) -> bool {
    lhs.len() < rhs.len() && lhs.iter().all(|element| rhs.binary_search(element).is_ok())
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, vector};

    use super::*;
    use crate::base::{Operation, Rotation};
    use crate::data::{Setting, operations_from_number};

    fn dihedral_three_operations() -> Operations {
        let r = matrix![
            0, -1, 0;
            1, -1, 0;
            0, 0, 1;
        ];
        let s = matrix![
            0, 1, 0;
            1, 0, 0;
            0, 0, 1;
        ];
        let rotations = [Rotation::identity(), r, r * r, s, r * s, r * r * s];
        rotations
            .into_iter()
            .map(|rotation| Operation::new(rotation, vector![0.0, 0.0, 0.0]))
            .collect()
    }

    #[test]
    fn test_subgroup_depths_without_ordering_assumption() {
        assert!(subgroup_depths(&[]).is_empty());

        let subgroups = vec![vec![0], vec![0, 1, 2, 3], vec![0, 2]];
        let depths = subgroup_depths(&subgroups);
        assert_eq!(depths[&[0, 1, 2, 3][..]], 0);
        assert_eq!(depths[&[0, 2][..]], 1);
        assert_eq!(depths[&[0][..]], 2);
    }

    #[test]
    fn test_enumerate_nonabelian_translationengleiche_subgroups() {
        let operations = dihedral_three_operations();
        let classes = enumerate_translationengleiche_subgroups(&operations, 1e-8).unwrap();

        assert_eq!(classes.len(), 4);
        assert_eq!(
            classes
                .iter()
                .map(|class| class.representative.operations.len())
                .collect::<Vec<_>>(),
            vec![6, 3, 2, 1]
        );
        let order_two = classes
            .iter()
            .find(|class| class.representative.operations.len() == 2)
            .unwrap();
        assert_eq!(order_two.conjugates.len(), 3);
        assert_eq!(order_two.representative.translationengleiche_index, 3);
        assert_eq!(order_two.representative.depth, 1);
        assert_eq!(order_two.normalizer_indices.len(), 2);
        assert_eq!(order_two.normalizer_permutations.len(), 2);
    }

    #[test]
    fn test_preserve_nonsymmorphic_subgroup_embeddings() {
        let operations = operations_from_number(19, Setting::Spglib, true).unwrap();
        let classes = enumerate_translationengleiche_subgroups(&operations, 1e-8).unwrap();

        assert_eq!(classes.len(), 5);
        assert_eq!(
            classes[0].representative.operation_indices,
            vec![0, 1, 2, 3]
        );
        assert_eq!(classes[0].representative.translationengleiche_index, 1);
        assert_eq!(classes[0].representative.depth, 0);
        assert_eq!(
            classes
                .last()
                .unwrap()
                .representative
                .translationengleiche_index,
            4
        );
        for class in &classes {
            for conjugate in &class.conjugates {
                for (operation, &parent_index) in conjugate
                    .subgroup
                    .operations
                    .iter()
                    .zip(&conjugate.subgroup.operation_indices)
                {
                    assert_eq!(operation.rotation, operations[parent_index].rotation);
                    assert_eq!(operation.translation, operations[parent_index].translation);
                }
            }
        }
    }

    #[test]
    fn test_reject_duplicate_rotations() {
        let operations = vec![Operation::identity(), Operation::identity()];
        assert_eq!(
            enumerate_translationengleiche_subgroups(&operations, 1e-8).unwrap_err(),
            MoyoError::InvalidPrimitiveOperationsError
        );
    }

    #[test]
    fn test_reject_operations_not_closed_modulo_translations() {
        let operations = vec![
            Operation::identity(),
            Operation::new(
                matrix![
                    -1, 0, 0;
                    0, 1, 0;
                    0, 0, 1;
                ],
                vector![0.0, 0.25, 0.0],
            ),
        ];
        assert_eq!(
            enumerate_translationengleiche_subgroups(&operations, 1e-8).unwrap_err(),
            MoyoError::InvalidPrimitiveOperationsError
        );
    }
}
