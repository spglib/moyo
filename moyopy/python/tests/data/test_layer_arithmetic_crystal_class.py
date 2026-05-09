from moyopy import LayerArithmeticCrystalClass


def test_layer_arithmetic_crystal_class_first():
    acc = LayerArithmeticCrystalClass(1)
    assert acc.arithmetic_number == 1
    assert acc.symbol == "p1"
    assert acc.geometric_crystal_class == "1"
    assert acc.bravais_class == "mp"
    assert acc.lattice_system == "Oblique"


def test_layer_arithmetic_crystal_class_centered():
    # arithmetic class 7 = c211 (centered rectangular)
    acc = LayerArithmeticCrystalClass(7)
    assert acc.symbol == "c211"
    assert acc.bravais_class == "oc"
    assert acc.lattice_system == "Rectangular"


def test_layer_arithmetic_crystal_class_last():
    acc = LayerArithmeticCrystalClass(43)
    assert acc.symbol == "p6/mmm"
    assert acc.geometric_crystal_class == "6/mmm"
    assert acc.bravais_class == "hp"
    assert acc.lattice_system == "Hexagonal"
