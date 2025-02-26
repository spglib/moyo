# ruff: noqa: E501
from moyopy import SpaceGroupType

map_hm_short_to_number = {
    # SpaceGroupType(number).hm_short is separated by a space like "F m -3 m"
    SpaceGroupType(number).hm_short.replace(" ", ""): number
    for number in range(1, 230 + 1)
}

print(map_hm_short_to_number.get("Fm-3m"))  # -> 225
print(map_hm_short_to_number.get("P6_3/m"))  # -> 176
print(map_hm_short_to_number.get("Fmmm"))  # -> 69
