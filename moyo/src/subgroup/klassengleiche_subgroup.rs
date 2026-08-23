use serde::Serialize;

use crate::base::{Linear, MoyoError, Operation, Operations, Transformation};
use crate::subgroup::affine_quotient::AffineQuotient;
use crate::subgroup::finite_group::FiniteGroup;

/// A Klassengleiche subgroup embedded in its parent space group.
///
/// A Klassengleiche subgroup has the same point group as its parent but a
/// finite-index translation sublattice. Its operations are given both in the
/// subgroup translation basis and in the parent primitive basis.
#[derive(Debug, Clone, Serialize)]
pub struct KlassengleicheSubgroup {
    /// Coset representatives in the subgroup primitive basis.
    pub operations: Operations,
    /// The same symmetry operations expressed in the parent primitive basis.
    pub parent_operations: Operations,
    /// Mapping from [`operations`](Self::operations) to the input parent operations.
    pub operation_indices: Vec<usize>,
    /// Basis transformation whose columns generate the translation sublattice.
    pub transformation: Linear,
    /// Klassengleiche index, equal to the translation-sublattice index.
    pub klassengleiche_index: usize,
}

/// A member of a parent-conjugacy class of embedded Klassengleiche subgroups.
#[derive(Debug, Clone, Serialize)]
pub struct KlassengleicheSubgroupConjugate {
    /// The conjugate subgroup embedded in the input parent.
    pub subgroup: KlassengleicheSubgroup,
    /// A parent-space-group operation whose inverse conjugates the
    /// representative to [`subgroup`](Self::subgroup).
    pub conjugator: Operation,
}

/// A parent-conjugacy class of embedded Klassengleiche subgroups.
#[derive(Debug, Clone, Serialize)]
pub struct KlassengleicheSubgroupConjugacyClass {
    /// Deterministically selected representative of this conjugacy class.
    pub representative: KlassengleicheSubgroup,
    /// Every distinct conjugate embedding, including the representative first.
    pub conjugates: Vec<KlassengleicheSubgroupConjugate>,
    /// Coset representatives of the normalizer modulo the subgroup translation
    /// lattice, expressed in the parent primitive basis.
    pub normalizer_operations: Operations,
    /// Permutations induced on the representative operations by the normalizer.
    ///
    /// This is aligned with [`normalizer_operations`](Self::normalizer_operations).
    /// For a permutation `p`, operation position `i` is mapped to `p[i]` by
    /// inverse conjugation with the corresponding normalizer operation.
    pub normalizer_permutations: Vec<Vec<usize>>,
}

/// Enumerate Klassengleiche subgroups with a prescribed translation sublattice.
///
/// `prim_operations` must contain exactly one representative for every point-
/// group element. `transformation` is expressed in the input primitive basis;
/// its columns generate the desired subgroup translation lattice. It must have
/// positive determinant and be preserved by every parent rotation.
///
/// The finite affine quotient `G/L` is formed explicitly. Its complements to
/// `T/L` correspond exactly to the Klassengleiche subgroups of `G` whose
/// translation lattice is `L`. The result contains all such subgroups,
/// partitioned into conjugacy classes under the parent space group.
pub fn enumerate_klassengleiche_subgroups(
    prim_operations: &Operations,
    transformation: &Linear,
    epsilon: f64,
) -> Result<Vec<KlassengleicheSubgroupConjugacyClass>, MoyoError> {
    FiniteGroup::from_operations(prim_operations, epsilon)?;
    let quotient = AffineQuotient::new(prim_operations, transformation, epsilon)?;
    let subgroups = quotient
        .group
        .enumerate_subgroups()
        .into_iter()
        .filter(|subgroup| quotient.is_complement(subgroup))
        .collect::<Vec<_>>();

    let mut classes = Vec::new();
    for conjugacy_class in quotient.group.conjugacy_classes(&subgroups) {
        let representative_indices = &conjugacy_class.representative;
        debug_assert!(
            conjugacy_class
                .conjugates
                .iter()
                .all(|conjugate| quotient.is_complement(&conjugate.subgroup))
        );
        let representative = quotient.make_subgroup(representative_indices);
        let conjugates = conjugacy_class
            .conjugates
            .into_iter()
            .map(|conjugate| KlassengleicheSubgroupConjugate {
                subgroup: quotient.make_subgroup(&conjugate.subgroup),
                conjugator: quotient.parent_operations[conjugate.conjugator_index].clone(),
            })
            .collect::<Vec<_>>();

        let normalizer_operations = conjugacy_class
            .normalizer_indices
            .into_iter()
            .map(|index| quotient.parent_operations[index].clone())
            .collect();
        classes.push(KlassengleicheSubgroupConjugacyClass {
            representative,
            conjugates,
            normalizer_operations,
            normalizer_permutations: conjugacy_class.normalizer_permutations,
        });
    }

    Ok(classes)
}

/// Enumerate Klassengleiche subgroups with a prescribed index.
///
/// Every index-`index` translation sublattice is generated once in Hermite
/// normal form. Sublattices not preserved by the parent point group cannot be
/// translation lattices of a Klassengleiche subgroup and are skipped. The
/// returned conjugacy classes retain each subgroup's exact embedding in the
/// input primitive basis.
pub fn enumerate_klassengleiche_subgroups_by_index(
    prim_operations: &Operations,
    index: usize,
    epsilon: f64,
) -> Result<Vec<KlassengleicheSubgroupConjugacyClass>, MoyoError> {
    FiniteGroup::from_operations(prim_operations, epsilon)?;

    let mut classes = Vec::new();
    for transformation in sublattice_transformations(index)? {
        if !is_sublattice_preserved(prim_operations, &transformation) {
            continue;
        }
        classes.append(&mut enumerate_klassengleiche_subgroups(
            prim_operations,
            &transformation,
            epsilon,
        )?);
    }
    Ok(classes)
}

fn is_sublattice_preserved(prim_operations: &Operations, transformation: &Linear) -> bool {
    let to_sublattice = Transformation::from_linear(*transformation);
    prim_operations
        .iter()
        .all(|operation| to_sublattice.transform_operation(operation).is_some())
}

fn sublattice_transformations(index: usize) -> Result<Vec<Linear>, MoyoError> {
    let index = i32::try_from(index).map_err(|_| MoyoError::InvalidSublatticeIndexError)?;
    if index < 1 {
        return Err(MoyoError::InvalidSublatticeIndexError);
    }

    let mut transformations = Vec::new();
    for a in 1..=index {
        if index % a != 0 {
            continue;
        }
        for c in 1..=index / a {
            if index % c != 0 || index % (a * c) != 0 {
                continue;
            }
            let f = index / (a * c);
            for b in 0..c {
                for d in 0..f {
                    for e in 0..f {
                        transformations.push(Linear::new(a, 0, 0, b, c, 0, d, e, f));
                    }
                }
            }
        }
    }
    Ok(transformations)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use nalgebra::{matrix, vector};

    use super::*;
    use crate::base::Rotation;
    use crate::data::{Setting, operations_from_number};
    use crate::subgroup::affine_quotient::translations_equivalent;

    #[test]
    fn test_enumerate_translation_sublattice() {
        let operations = vec![Operation::identity()];
        let transformation = matrix![
            2, 0, 0;
            0, 1, 0;
            0, 0, 1;
        ];
        let classes =
            enumerate_klassengleiche_subgroups(&operations, &transformation, 1e-8).unwrap();

        assert_eq!(classes.len(), 1);
        let subgroup = &classes[0].representative;
        assert_eq!(subgroup.operations.len(), 1);
        assert_eq!(subgroup.operation_indices, vec![0]);
        assert_eq!(subgroup.transformation, transformation);
        assert_eq!(subgroup.klassengleiche_index, 2);
        assert_eq!(classes[0].conjugates.len(), 1);
        assert_eq!(classes[0].normalizer_operations.len(), 2);
        assert_eq!(classes[0].normalizer_permutations, vec![vec![0], vec![0]]);
    }

    #[test]
    fn test_enumerate_klassengleiche_subgroups_by_index() {
        let operations = vec![Operation::identity()];
        let classes = enumerate_klassengleiche_subgroups_by_index(&operations, 2, 1e-8).unwrap();

        assert_eq!(classes.len(), 7);
        assert!(classes.iter().all(|class| {
            class.representative.klassengleiche_index == 2
                && class.representative.operations.len() == 1
        }));
        assert_eq!(
            classes
                .iter()
                .map(|class| class.representative.transformation.as_slice().to_vec())
                .collect::<BTreeSet<_>>()
                .len(),
            7
        );
    }

    #[test]
    fn test_reject_invalid_sublattice_index() {
        let operations = vec![Operation::identity()];
        for index in [0, usize::MAX] {
            assert_eq!(
                enumerate_klassengleiche_subgroups_by_index(&operations, index, 1e-8).unwrap_err(),
                MoyoError::InvalidSublatticeIndexError
            );
        }
    }

    #[test]
    fn test_enumerate_distinct_inversion_complements() {
        let operations = vec![
            Operation::identity(),
            Operation::new(-Rotation::identity(), vector![0.0, 0.0, 0.0]),
        ];
        let transformation = matrix![
            2, 0, 0;
            0, 1, 0;
            0, 0, 1;
        ];
        let classes =
            enumerate_klassengleiche_subgroups(&operations, &transformation, 1e-8).unwrap();

        assert_eq!(classes.len(), 2);
        for class in &classes {
            assert_eq!(class.representative.operations.len(), 2);
            assert_eq!(class.representative.operation_indices, vec![0, 1]);
            assert_eq!(class.representative.klassengleiche_index, 2);
            assert_eq!(class.conjugates.len(), 1);
            assert_eq!(class.normalizer_operations.len(), 4);
        }
        let inversion_translations = classes
            .iter()
            .map(|class| class.representative.operations[1].translation)
            .collect::<Vec<_>>();
        assert!(inversion_translations.contains(&vector![0.0, 0.0, 0.0]));
        assert!(inversion_translations.contains(&vector![0.5, 0.0, 0.0]));
        let parent_inversion_translations = classes
            .iter()
            .map(|class| class.representative.parent_operations[1].translation)
            .collect::<Vec<_>>();
        assert!(parent_inversion_translations.contains(&vector![0.0, 0.0, 0.0]));
        assert!(parent_inversion_translations.contains(&vector![1.0, 0.0, 0.0]));
    }

    #[test]
    fn test_partition_nontrivial_parent_conjugacy() {
        let fourfold = matrix![
            0, -1, 0;
            1, 0, 0;
            0, 0, 1;
        ];
        let rotations = [
            Rotation::identity(),
            fourfold,
            fourfold * fourfold,
            fourfold * fourfold * fourfold,
        ];
        let operations = rotations
            .into_iter()
            .map(|rotation| Operation::new(rotation, vector![0.0, 0.0, 0.0]))
            .collect();
        let transformation = matrix![
            2, 0, 0;
            0, 2, 0;
            0, 0, 1;
        ];
        let classes =
            enumerate_klassengleiche_subgroups(&operations, &transformation, 1e-8).unwrap();

        assert_eq!(classes.len(), 2);
        for class in &classes {
            assert_eq!(class.conjugates.len(), 2);
            assert_eq!(class.normalizer_operations.len(), 8);
            assert_eq!(class.representative.operation_indices, vec![0, 1, 2, 3]);
            assert_eq!(class.normalizer_permutations.len(), 8);
        }
    }

    #[test]
    fn test_preserve_nonsymmorphic_klassengleiche_embeddings() {
        let operations = operations_from_number(19, Setting::Spglib, true).unwrap();
        let transformation = matrix![
            3, 0, 0;
            0, 1, 0;
            0, 0, 1;
        ];
        let classes =
            enumerate_klassengleiche_subgroups(&operations, &transformation, 1e-8).unwrap();

        assert!(!classes.is_empty());
        for class in &classes {
            for conjugate in &class.conjugates {
                assert_eq!(conjugate.subgroup.operation_indices, vec![0, 1, 2, 3]);
                for ((operation, parent_operation), &parent_index) in conjugate
                    .subgroup
                    .operations
                    .iter()
                    .zip(&conjugate.subgroup.parent_operations)
                    .zip(&conjugate.subgroup.operation_indices)
                {
                    assert_eq!(parent_operation.rotation, operations[parent_index].rotation);
                    let transformed = Transformation::from_linear(transformation)
                        .transform_operation(parent_operation)
                        .unwrap();
                    assert_eq!(operation.rotation, transformed.rotation);
                    assert!(translations_equivalent(
                        &operation.translation,
                        &transformed.translation,
                        1e-8
                    ));
                }
            }
        }
    }

    #[test]
    fn test_reject_sublattice_not_preserved_by_point_group() {
        let fourfold = matrix![
            0, -1, 0;
            1, 0, 0;
            0, 0, 1;
        ];
        let rotations = [
            Rotation::identity(),
            fourfold,
            fourfold * fourfold,
            fourfold * fourfold * fourfold,
        ];
        let operations = rotations
            .into_iter()
            .map(|rotation| Operation::new(rotation, vector![0.0, 0.0, 0.0]))
            .collect();
        let transformation = matrix![
            2, 0, 0;
            0, 1, 0;
            0, 0, 1;
        ];

        assert_eq!(
            enumerate_klassengleiche_subgroups(&operations, &transformation, 1e-8).unwrap_err(),
            MoyoError::InvalidSublatticeTransformationError
        );
    }

    #[test]
    fn test_reject_nonpositive_sublattice_determinant() {
        let operations = vec![Operation::identity()];
        for transformation in [
            matrix![
                0, 0, 0;
                0, 1, 0;
                0, 0, 1;
            ],
            matrix![
                -1, 0, 0;
                0, 1, 0;
                0, 0, 1;
            ],
            matrix![
                i32::MAX, 0, 0;
                0, 2, 0;
                0, 0, 1;
            ],
        ] {
            assert_eq!(
                enumerate_klassengleiche_subgroups(&operations, &transformation, 1e-8).unwrap_err(),
                MoyoError::InvalidSublatticeTransformationError
            );
        }
    }
}
