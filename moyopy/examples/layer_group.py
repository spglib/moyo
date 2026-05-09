from math import sqrt

import moyopy

# MoS2 1H monolayer (layer group 78, p-6m2).
# In-plane lattice constant a ~ 3.16 A; intra-layer S-S spacing ~ 3.12 A.
# The third basis vector is the aperiodic stacking direction; its length is
# arbitrary (moyo does not rescale c). See moyopy/docs/layer_standardization.md.
a = 3.16
c = 20.0
basis = [
    [a, 0.0, 0.0],
    [-a / 2.0, a * sqrt(3.0) / 2.0, 0.0],
    [0.0, 0.0, c],
]

# Place the Mo layer at z = 0.5 so the two S atoms straddle it within [0, 1).
# In 1H MoS2 the two S atoms sit directly above/below the same in-plane site,
# giving trigonal-prismatic Mo coordination (no inversion through Mo).
zS = 1.56 / c
positions = [
    [0.0, 0.0, 0.5],  # Mo
    [1.0 / 3.0, 2.0 / 3.0, 0.5 + zS],  # S top
    [1.0 / 3.0, 2.0 / 3.0, 0.5 - zS],  # S bottom
]
numbers = [42, 16, 16]
cell = moyopy.Cell(basis, positions, numbers)

dataset = moyopy.MoyoLayerDataset(cell, symprec=1e-4)
assert dataset.number == 78
assert dataset.hall_number == 114
assert dataset.pearson_symbol == "hp3"
assert dataset.orbits == [0, 1, 1]
assert dataset.wyckoffs == ["a", "e", "e"]
assert dataset.site_symmetry_symbols == ["-6m2", "3m.", "3m."]
# |c| is preserved and c lies on Cartesian z in the standardized cell.
std_basis = dataset.std_cell.basis
assert abs(std_basis[2][0]) < 1e-6
assert abs(std_basis[2][1]) < 1e-6
assert abs(std_basis[2][2] - c) < 1e-6

# MoyoLayerDataset round-trips through JSON.
dataset2 = moyopy.MoyoLayerDataset.deserialize_json(dataset.serialize_json())
assert dataset2.number == dataset.number
assert dataset2.hall_number == dataset.hall_number

print(dataset)
