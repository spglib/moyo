use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};

use serde::Serialize;

use crate::base::{MoyoError, Operations, Rotation};

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

    let mut remaining = subgroups.iter().cloned().collect::<BTreeSet<_>>();
    let mut classes = Vec::new();
    for representative_indices in &subgroups {
        if !remaining.remove(representative_indices) {
            continue;
        }

        let mut conjugate_witnesses =
            BTreeMap::from([(representative_indices.clone(), finite_group.identity)]);
        for conjugator in 0..finite_group.order() {
            let conjugate = finite_group.conjugate(representative_indices, conjugator);
            conjugate_witnesses.entry(conjugate).or_insert(conjugator);
        }
        for conjugate in conjugate_witnesses.keys() {
            remaining.remove(conjugate);
        }

        let representative = make_translationengleiche_subgroup(
            representative_indices,
            prim_operations,
            finite_group.order(),
            &depths,
        );
        let identity_index = finite_group.identity;
        let mut conjugates = conjugate_witnesses
            .into_iter()
            .map(
                |(indices, conjugator_index)| TranslationengleicheSubgroupConjugate {
                    subgroup: make_translationengleiche_subgroup(
                        &indices,
                        prim_operations,
                        finite_group.order(),
                        &depths,
                    ),
                    conjugator_index,
                },
            )
            .collect::<Vec<_>>();
        conjugates.sort_by_key(|conjugate| {
            (
                conjugate.subgroup.operation_indices != *representative_indices,
                conjugate.subgroup.operation_indices.clone(),
            )
        });
        debug_assert_eq!(conjugates[0].conjugator_index, identity_index);

        let (normalizer_indices, normalizer_permutations) =
            finite_group.normalizer_action(representative_indices);
        classes.push(TranslationengleicheSubgroupConjugacyClass {
            representative,
            conjugates,
            normalizer_indices,
            normalizer_permutations,
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

struct FiniteGroup {
    table: Vec<Vec<usize>>,
    inverses: Vec<usize>,
    identity: usize,
}

impl FiniteGroup {
    fn from_operations(operations: &Operations, epsilon: f64) -> Result<Self, MoyoError> {
        if operations.is_empty() || !epsilon.is_finite() || epsilon <= 0.0 {
            return Err(MoyoError::InvalidPrimitiveOperationsError);
        }

        let mut rotation_indices = HashMap::with_capacity(operations.len());
        for (index, operation) in operations.iter().enumerate() {
            if rotation_indices.insert(operation.rotation, index).is_some() {
                return Err(MoyoError::InvalidPrimitiveOperationsError);
            }
        }

        let identity = *rotation_indices
            .get(&Rotation::identity())
            .ok_or(MoyoError::InvalidPrimitiveOperationsError)?;
        let mut table = vec![vec![0; operations.len()]; operations.len()];
        for (lhs_index, lhs) in operations.iter().enumerate() {
            for (rhs_index, rhs) in operations.iter().enumerate() {
                let product = lhs.clone() * rhs.clone();
                let product_index = *rotation_indices
                    .get(&product.rotation)
                    .ok_or(MoyoError::InvalidPrimitiveOperationsError)?;
                let translation_difference =
                    product.translation - operations[product_index].translation;
                if translation_difference
                    .iter()
                    .any(|&value| (value - value.round()).abs() >= epsilon)
                {
                    return Err(MoyoError::InvalidPrimitiveOperationsError);
                }
                table[lhs_index][rhs_index] = product_index;
            }
        }

        let inverses = (0..operations.len())
            .map(|element| {
                (0..operations.len())
                    .find(|&candidate| {
                        table[element][candidate] == identity
                            && table[candidate][element] == identity
                    })
                    .ok_or(MoyoError::InvalidPrimitiveOperationsError)
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self {
            table,
            inverses,
            identity,
        })
    }

    fn order(&self) -> usize {
        self.table.len()
    }

    fn enumerate_subgroups(&self) -> Vec<Vec<usize>> {
        let trivial = vec![self.identity];
        let mut seen = BTreeSet::from([trivial.clone()]);
        let mut queue = VecDeque::from([trivial]);

        while let Some(subgroup) = queue.pop_front() {
            for generator in 0..self.order() {
                if subgroup.binary_search(&generator).is_ok() {
                    continue;
                }
                let generated = self.generated_subgroup(&subgroup, generator);
                if seen.insert(generated.clone()) {
                    queue.push_back(generated);
                }
            }
        }

        let mut subgroups = seen.into_iter().collect::<Vec<_>>();
        subgroups.sort_by(|lhs, rhs| rhs.len().cmp(&lhs.len()).then_with(|| lhs.cmp(rhs)));
        subgroups
    }

    fn generated_subgroup(&self, subgroup: &[usize], generator: usize) -> Vec<usize> {
        let mut members = vec![false; self.order()];
        members[self.identity] = true;
        members[generator] = true;
        for &element in subgroup {
            members[element] = true;
        }

        loop {
            let elements = members
                .iter()
                .enumerate()
                .filter_map(|(index, &present)| present.then_some(index))
                .collect::<Vec<_>>();
            let mut changed = false;
            for &lhs in &elements {
                for &rhs in &elements {
                    let product = self.table[lhs][rhs];
                    if !members[product] {
                        members[product] = true;
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }

        members
            .into_iter()
            .enumerate()
            .filter_map(|(index, present)| present.then_some(index))
            .collect()
    }

    fn conjugate(&self, subgroup: &[usize], conjugator: usize) -> Vec<usize> {
        let inverse = self.inverses[conjugator];
        let mut conjugate = subgroup
            .iter()
            .map(|&element| self.table[self.table[inverse][element]][conjugator])
            .collect::<Vec<_>>();
        conjugate.sort_unstable();
        conjugate
    }

    fn normalizer_action(&self, subgroup: &[usize]) -> (Vec<usize>, Vec<Vec<usize>>) {
        let local_indices = subgroup
            .iter()
            .enumerate()
            .map(|(local, &parent)| (parent, local))
            .collect::<HashMap<_, _>>();
        let mut normalizer_indices = Vec::new();
        let mut normalizer_permutations = Vec::new();

        for conjugator in 0..self.order() {
            if self.conjugate(subgroup, conjugator) != subgroup {
                continue;
            }
            let inverse = self.inverses[conjugator];
            let permutation = subgroup
                .iter()
                .map(|&element| {
                    let image = self.table[self.table[inverse][element]][conjugator];
                    local_indices[&image]
                })
                .collect();
            normalizer_indices.push(conjugator);
            normalizer_permutations.push(permutation);
        }

        (normalizer_indices, normalizer_permutations)
    }
}

fn subgroup_depths(subgroups: &[Vec<usize>]) -> BTreeMap<Vec<usize>, usize> {
    let mut depths = BTreeMap::new();
    depths.insert(subgroups[0].clone(), 0);

    for subgroup in subgroups.iter().skip(1) {
        let covering_supergroups = subgroups.iter().filter(|candidate| {
            is_proper_subset(subgroup, candidate)
                && !subgroups.iter().any(|between| {
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
    use crate::base::Operation;
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
    fn test_lift_nonsymmorphic_translationengleiche_subgroups() {
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
