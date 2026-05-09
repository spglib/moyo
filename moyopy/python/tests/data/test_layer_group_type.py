from moyopy import LayerGroupType


def test_layer_group_type_lg1():
    lgt = LayerGroupType(1)
    assert lgt.number == 1
    assert lgt.hm_short == "p 1"
    assert lgt.hm_full == "p 1"
    assert lgt.arithmetic_number == 1
    assert lgt.arithmetic_symbol == "p1"
    assert lgt.geometric_crystal_class == "1"
    assert lgt.bravais_class == "mp"
    assert lgt.lattice_system == "Oblique"


def test_layer_group_type_lg10_centered():
    # LG 10 = c211, the smallest centered (oc) layer group.
    lgt = LayerGroupType(10)
    assert lgt.number == 10
    assert lgt.arithmetic_symbol == "c211"
    assert lgt.geometric_crystal_class == "2"
    assert lgt.bravais_class == "oc"
    assert lgt.lattice_system == "Rectangular"


def test_layer_group_type_lg80():
    lgt = LayerGroupType(80)
    assert lgt.number == 80
    assert lgt.arithmetic_number == 43
    assert lgt.arithmetic_symbol == "p6/mmm"
    assert lgt.geometric_crystal_class == "6/mmm"
    assert lgt.bravais_class == "hp"
    assert lgt.lattice_system == "Hexagonal"
