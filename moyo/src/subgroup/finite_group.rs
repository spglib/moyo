use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};

use crate::base::{MoyoError, Operations, Rotation};

pub(super) struct FiniteGroup {
    table: Vec<Vec<usize>>,
    inverses: Vec<usize>,
    identity: usize,
}

pub(super) struct FiniteSubgroupConjugate {
    pub(super) subgroup: Vec<usize>,
    pub(super) conjugator_index: usize,
}

pub(super) struct FiniteSubgroupConjugacyClass {
    pub(super) representative: Vec<usize>,
    pub(super) conjugates: Vec<FiniteSubgroupConjugate>,
    pub(super) normalizer_indices: Vec<usize>,
    pub(super) normalizer_permutations: Vec<Vec<usize>>,
}

impl FiniteGroup {
    pub(super) fn from_table(table: Vec<Vec<usize>>) -> Option<Self> {
        let order = table.len();
        if order == 0
            || table
                .iter()
                .any(|row| row.len() != order || row.iter().any(|&element| element >= order))
        {
            return None;
        }

        let identity = (0..order).find(|&candidate| {
            (0..order).all(|element| {
                table[candidate][element] == element && table[element][candidate] == element
            })
        })?;
        let inverses = (0..order)
            .map(|element| {
                (0..order).find(|&candidate| {
                    table[element][candidate] == identity && table[candidate][element] == identity
                })
            })
            .collect::<Option<Vec<_>>>()?;

        Some(Self {
            table,
            inverses,
            identity,
        })
    }

    pub(super) fn from_operations(
        operations: &Operations,
        epsilon: f64,
    ) -> Result<Self, MoyoError> {
        if operations.is_empty() || !epsilon.is_finite() || epsilon <= 0.0 {
            return Err(MoyoError::InvalidPrimitiveOperationsError);
        }

        let mut rotation_indices = HashMap::with_capacity(operations.len());
        for (index, operation) in operations.iter().enumerate() {
            if rotation_indices.insert(operation.rotation, index).is_some() {
                return Err(MoyoError::InvalidPrimitiveOperationsError);
            }
        }

        if !rotation_indices.contains_key(&Rotation::identity()) {
            return Err(MoyoError::InvalidPrimitiveOperationsError);
        }
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

        Self::from_table(table).ok_or(MoyoError::InvalidPrimitiveOperationsError)
    }

    pub(super) fn identity(&self) -> usize {
        self.identity
    }

    pub(super) fn order(&self) -> usize {
        self.table.len()
    }

    pub(super) fn enumerate_subgroups(&self) -> Vec<Vec<usize>> {
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

    pub(super) fn conjugate(&self, subgroup: &[usize], conjugator: usize) -> Vec<usize> {
        let inverse = self.inverses[conjugator];
        let mut conjugate = subgroup
            .iter()
            .map(|&element| self.table[self.table[inverse][element]][conjugator])
            .collect::<Vec<_>>();
        conjugate.sort_unstable();
        conjugate
    }

    pub(super) fn conjugacy_classes(
        &self,
        subgroups: &[Vec<usize>],
    ) -> Vec<FiniteSubgroupConjugacyClass> {
        let mut remaining = subgroups.iter().cloned().collect::<BTreeSet<_>>();
        let mut classes = Vec::new();
        for representative in subgroups {
            if !remaining.remove(representative) {
                continue;
            }

            let mut conjugate_witnesses =
                BTreeMap::from([(representative.clone(), self.identity())]);
            for conjugator in 0..self.order() {
                let conjugate = self.conjugate(representative, conjugator);
                conjugate_witnesses.entry(conjugate).or_insert(conjugator);
            }
            for conjugate in conjugate_witnesses.keys() {
                remaining.remove(conjugate);
            }

            let mut conjugates = conjugate_witnesses
                .into_iter()
                .map(|(subgroup, conjugator_index)| FiniteSubgroupConjugate {
                    subgroup,
                    conjugator_index,
                })
                .collect::<Vec<_>>();
            conjugates.sort_by_key(|conjugate| {
                (
                    conjugate.subgroup != *representative,
                    conjugate.subgroup.clone(),
                )
            });
            debug_assert_eq!(conjugates[0].conjugator_index, self.identity());

            let (normalizer_indices, normalizer_permutations) =
                self.normalizer_action(representative);
            classes.push(FiniteSubgroupConjugacyClass {
                representative: representative.clone(),
                conjugates,
                normalizer_indices,
                normalizer_permutations,
            });
        }
        classes
    }

    pub(super) fn normalizer_action(&self, subgroup: &[usize]) -> (Vec<usize>, Vec<Vec<usize>>) {
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
