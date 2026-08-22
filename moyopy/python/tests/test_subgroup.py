import pytest

from moyopy import (
    KlassengleicheSubgroup,
    KlassengleicheSubgroupConjugacyClass,
    KlassengleicheSubgroupConjugate,
    TranslationengleicheSubgroup,
    TranslationengleicheSubgroupConjugacyClass,
    TranslationengleicheSubgroupConjugate,
    enumerate_klassengleiche_subgroups,
    enumerate_klassengleiche_subgroups_by_index,
    enumerate_translationengleiche_subgroups,
    operations_from_number,
)


def test_enumerate_klassengleiche_subgroups():
    operations = operations_from_number(1, primitive=True)  # P1
    transformation = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    classes = enumerate_klassengleiche_subgroups(
        operations.rotations, operations.translations, transformation
    )

    assert len(classes) == 1
    conjugacy_class = classes[0]
    assert isinstance(conjugacy_class, KlassengleicheSubgroupConjugacyClass)
    subgroup = conjugacy_class.representative
    assert isinstance(subgroup, KlassengleicheSubgroup)
    assert len(subgroup.operations) == 1
    assert len(subgroup.parent_operations) == 1
    assert subgroup.operation_indices == [0]
    assert subgroup.transformation == transformation
    assert subgroup.klassengleiche_index == 2
    assert len(conjugacy_class.conjugates) == 1
    assert isinstance(conjugacy_class.conjugates[0], KlassengleicheSubgroupConjugate)
    assert len(conjugacy_class.conjugates[0].conjugator) == 1
    assert len(conjugacy_class.normalizer_operations) == 2
    assert conjugacy_class.normalizer_permutations == [[0], [0]]


def test_enumerate_klassengleiche_subgroups_by_index():
    operations = operations_from_number(1, primitive=True)  # P1
    classes = enumerate_klassengleiche_subgroups_by_index(
        operations.rotations, operations.translations, 2
    )

    assert len(classes) == 7
    assert all(item.representative.klassengleiche_index == 2 for item in classes)
    assert len({str(item.representative.transformation) for item in classes}) == 7


def test_enumerate_translationengleiche_subgroups():
    operations = operations_from_number(149, primitive=True)  # P312, point group 32
    classes = enumerate_translationengleiche_subgroups(
        operations.rotations, operations.translations
    )

    assert len(classes) == 4
    assert all(isinstance(item, TranslationengleicheSubgroupConjugacyClass) for item in classes)
    assert [len(item.representative.operations) for item in classes] == [6, 3, 2, 1]
    order_two = next(item for item in classes if len(item.representative.operations) == 2)
    assert isinstance(order_two.representative, TranslationengleicheSubgroup)
    assert order_two.representative.translationengleiche_index == 3
    assert order_two.representative.depth == 1
    assert len(order_two.conjugates) == 3
    assert all(
        isinstance(item, TranslationengleicheSubgroupConjugate) for item in order_two.conjugates
    )
    assert len(order_two.normalizer_indices) == 2
    assert len(order_two.normalizer_permutations) == 2


def test_enumerate_translationengleiche_subgroups_rejects_mismatched_operations():
    with pytest.raises(ValueError, match="must have the same length"):
        enumerate_translationengleiche_subgroups([[[1, 0, 0], [0, 1, 0], [0, 0, 1]]], [])


def test_enumerate_klassengleiche_subgroups_rejects_mismatched_operations():
    with pytest.raises(ValueError, match="must have the same length"):
        enumerate_klassengleiche_subgroups(
            [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
            [],
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        )
