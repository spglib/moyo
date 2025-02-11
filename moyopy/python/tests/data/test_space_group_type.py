from moyopy import SpaceGroupType


def test_space_group_type():
    s1 = SpaceGroupType(221)
    assert s1.number == 221
    assert s1.hm_short == "P m -3 m"
    assert s1.arithmetic_number == 71
    assert s1.arithmetic_symbol == "m-3mP"
    assert s1.geometric_crystal_class == "m-3m"
    assert s1.crystal_system == "Cubic"
    assert s1.bravais_class == "cP"
    assert s1.lattice_system == "Cubic"
    assert s1.crystal_family == "Cubic"

    s2 = SpaceGroupType(167)
    assert s2.number == 167
    assert s2.hm_short == "R -3 c"
    assert s2.arithmetic_number == 50
    assert s2.arithmetic_symbol == "-3mR"
    assert s2.geometric_crystal_class == "-3m"
    assert s2.crystal_system == "Trigonal"
    assert s2.bravais_class == "hR"
    assert s2.lattice_system == "Rhombohedral"
    assert s2.crystal_family == "Hexagonal"
