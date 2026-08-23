use std::collections::{BTreeMap, BTreeSet, HashMap};

use itertools::iproduct;
use nalgebra::Vector3;
use serde::Serialize;

use crate::base::{Linear, MoyoError, Operation, Operations, Rotation, Transformation};
use crate::math::SNF;
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

    let mut remaining = subgroups.iter().cloned().collect::<BTreeSet<_>>();
    let mut classes = Vec::new();
    for representative_indices in &subgroups {
        if !remaining.remove(representative_indices) {
            continue;
        }

        let mut conjugate_witnesses =
            BTreeMap::from([(representative_indices.clone(), quotient.group.identity())]);
        for conjugator in 0..quotient.group.order() {
            let conjugate = quotient.group.conjugate(representative_indices, conjugator);
            debug_assert!(quotient.is_complement(&conjugate));
            conjugate_witnesses.entry(conjugate).or_insert(conjugator);
        }
        for conjugate in conjugate_witnesses.keys() {
            remaining.remove(conjugate);
        }

        let representative = quotient.make_subgroup(representative_indices);
        let identity_index = quotient.group.identity();
        let mut conjugate_witnesses = conjugate_witnesses.into_iter().collect::<Vec<_>>();
        conjugate_witnesses
            .sort_by_key(|(indices, _)| (indices != representative_indices, indices.clone()));
        debug_assert_eq!(conjugate_witnesses[0].1, identity_index);
        let conjugates = conjugate_witnesses
            .into_iter()
            .map(
                |(indices, conjugator_index)| KlassengleicheSubgroupConjugate {
                    subgroup: quotient.make_subgroup(&indices),
                    conjugator: quotient.parent_operations[conjugator_index].clone(),
                },
            )
            .collect::<Vec<_>>();

        let (normalizer_indices, normalizer_permutations) =
            quotient.group.normalizer_action(representative_indices);
        let normalizer_operations = normalizer_indices
            .into_iter()
            .map(|index| quotient.parent_operations[index].clone())
            .collect();
        classes.push(KlassengleicheSubgroupConjugacyClass {
            representative,
            conjugates,
            normalizer_operations,
            normalizer_permutations,
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

/// Finite affine quotient `G/L` used to enumerate Klassengleiche subgroups.
///
/// `G` is the parent space group generated by the input primitive operations
/// and integer translations, and `L` is the point-group-invariant translation
/// sublattice generated by the columns of `transformation`. The quotient has
/// order `|G/T| * [T:L]`, where `T` is the parent translation subgroup.
///
/// `operations` and `parent_operations` are aligned representatives of the
/// same quotient elements. The former are expressed in the `L` basis with
/// translations reduced modulo integers; the latter are expressed in the
/// parent primitive basis. `point_indices` maps each quotient element to its
/// input point-operation representative.
///
/// A complement to `T/L` in this quotient has an inverse image in `G` whose
/// translation subgroup is exactly `L`, and therefore determines a
/// Klassengleiche subgroup with the prescribed translation lattice.
struct AffineQuotient {
    group: FiniteGroup,
    operations: Operations,
    parent_operations: Operations,
    point_indices: Vec<usize>,
    transformation: Linear,
    sublattice_index: usize,
    point_order: usize,
}

impl AffineQuotient {
    fn new(
        prim_operations: &Operations,
        transformation: &Linear,
        epsilon: f64,
    ) -> Result<Self, MoyoError> {
        let determinant = exact_determinant(transformation);
        if determinant <= 0 || determinant > i32::MAX as i128 {
            return Err(MoyoError::InvalidSublatticeTransformationError);
        }
        let sublattice_index = determinant as usize;
        let lattice_points = lattice_points(transformation);
        if lattice_points.len() != sublattice_index {
            return Err(MoyoError::InvalidSublatticeTransformationError);
        }

        let to_sublattice = Transformation::from_linear(*transformation);
        let quotient_order = prim_operations.len() * sublattice_index;
        let mut operations = Vec::with_capacity(quotient_order);
        let mut parent_operations = Vec::with_capacity(quotient_order);
        let mut point_indices = Vec::with_capacity(quotient_order);
        for (point_index, operation) in prim_operations.iter().enumerate() {
            for lattice_point in &lattice_points {
                let parent_operation = Operation::new(
                    operation.rotation,
                    operation.translation + lattice_point.map(|element| element as f64),
                );
                let mut transformed = to_sublattice
                    .transform_operation(&parent_operation)
                    .ok_or(MoyoError::InvalidSublatticeTransformationError)?;
                transformed.translation = transformed
                    .translation
                    .map(|element| reduce_mod_one(element, epsilon));
                operations.push(transformed);
                parent_operations.push(parent_operation);
                point_indices.push(point_index);
            }
        }

        let table = affine_cayley_table(&operations, epsilon)
            .ok_or(MoyoError::InvalidSublatticeTransformationError)?;
        let group = FiniteGroup::from_table(table)
            .ok_or(MoyoError::InvalidSublatticeTransformationError)?;

        Ok(Self {
            group,
            operations,
            parent_operations,
            point_indices,
            transformation: *transformation,
            sublattice_index,
            point_order: prim_operations.len(),
        })
    }

    fn is_complement(&self, subgroup: &[usize]) -> bool {
        if subgroup.len() != self.point_order {
            return false;
        }

        let mut covered = vec![false; self.point_order];
        for &index in subgroup {
            let point_index = self.point_indices[index];
            if covered[point_index] {
                return false;
            }
            covered[point_index] = true;
        }
        covered.into_iter().all(|value| value)
    }

    fn make_subgroup(&self, indices: &[usize]) -> KlassengleicheSubgroup {
        KlassengleicheSubgroup {
            operations: indices
                .iter()
                .map(|&index| self.operations[index].clone())
                .collect(),
            parent_operations: indices
                .iter()
                .map(|&index| self.parent_operations[index].clone())
                .collect(),
            operation_indices: indices
                .iter()
                .map(|&index| self.point_indices[index])
                .collect(),
            transformation: self.transformation,
            klassengleiche_index: self.sublattice_index,
        }
    }
}

fn lattice_points(transformation: &Linear) -> Vec<Vector3<i32>> {
    let snf = SNF::new(transformation);
    let linear_inverse = snf
        .l
        .map(|element| element as f64)
        .try_inverse()
        .expect("Smith transformation is unimodular")
        .map(|element| element.round() as i32);

    iproduct!(0..snf.d[(0, 0)], 0..snf.d[(1, 1)], 0..snf.d[(2, 2)])
        .map(|(f0, f1, f2)| linear_inverse * Vector3::new(f0, f1, f2))
        .collect()
}

fn exact_determinant(matrix: &Linear) -> i128 {
    let value = |row, column| matrix[(row, column)] as i128;
    value(0, 0) * (value(1, 1) * value(2, 2) - value(1, 2) * value(2, 1))
        - value(0, 1) * (value(1, 0) * value(2, 2) - value(1, 2) * value(2, 0))
        + value(0, 2) * (value(1, 0) * value(2, 1) - value(1, 1) * value(2, 0))
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

fn affine_cayley_table(operations: &Operations, epsilon: f64) -> Option<Vec<Vec<usize>>> {
    let mut rotation_indices = HashMap::<Rotation, Vec<usize>>::new();
    for (index, operation) in operations.iter().enumerate() {
        let candidates = rotation_indices.entry(operation.rotation).or_default();
        if candidates.iter().any(|&candidate| {
            translations_equivalent(
                &operation.translation,
                &operations[candidate].translation,
                epsilon,
            )
        }) {
            return None;
        }
        candidates.push(index);
    }

    let mut table = vec![vec![0; operations.len()]; operations.len()];
    for (lhs_index, lhs) in operations.iter().enumerate() {
        for (rhs_index, rhs) in operations.iter().enumerate() {
            let product = lhs.clone() * rhs.clone();
            let product_index =
                rotation_indices
                    .get(&product.rotation)?
                    .iter()
                    .find(|&&candidate| {
                        translations_equivalent(
                            &product.translation,
                            &operations[candidate].translation,
                            epsilon,
                        )
                    })?;
            table[lhs_index][rhs_index] = *product_index;
        }
    }
    Some(table)
}

fn translations_equivalent(lhs: &Vector3<f64>, rhs: &Vector3<f64>, epsilon: f64) -> bool {
    (lhs - rhs)
        .iter()
        .all(|&value| (value - value.round()).abs() < epsilon)
}

fn reduce_mod_one(value: f64, epsilon: f64) -> f64 {
    let reduced = value.rem_euclid(1.0);
    if reduced < epsilon || 1.0 - reduced < epsilon {
        0.0
    } else {
        reduced
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, vector};

    use super::*;
    use crate::data::{Setting, operations_from_number};

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
