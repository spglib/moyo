# ruff: noqa: E501
from moyopy import SpaceGroupType

s = SpaceGroupType(15)  # ITA space group number (1 - 230)
assert s.hm_short == "C 2/c"
assert s.crystal_system == "Monoclinic"

print(s)
# -> PySpaceGroupType { number: 15, hm_short: "C 2/c", hm_full: "C 1 2/c 1", arithmetic_number: 8, arithmetic_symbol: "2/mC", geometric_crystal_class: "2/m", crystal_system: "Monoclinic", bravais_class: "mC", lattice_system: "Monoclinic", crystal_family: "Monoclinic" }
