from moyopy import MagneticSpaceGroupType


def test_magnetic_space_group_type():
    msgt = MagneticSpaceGroupType(1262)
    assert msgt.bns_number == "151.32"
    assert msgt.og_number == "153.4.1270"
    assert msgt.construct_type == 4
