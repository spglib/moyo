import pytest

from moyopy import (
    TranslationengleicheSubgroup,
    TranslationengleicheSubgroupConjugacyClass,
    TranslationengleicheSubgroupConjugate,
    enumerate_translationengleiche_subgroups,
    operations_from_number,
)


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
