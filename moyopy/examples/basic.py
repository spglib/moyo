from math import sqrt

import moyopy

# https://next-gen.materialsproject.org/materials/mp-560588
a = 3.81
c = 6.24
basis = [
    [a, 0.0, 0.0],
    [-a / 2.0, a * sqrt(3.0) / 2.0, 0.0],
    [0.0, 0.0, c],
]
z1_2b = 0.00014
z2_2b = 0.37486
positions = [
    # 2b
    [1 / 3, 2 / 3, z1_2b],
    [2 / 3, 1 / 3, z1_2b + 0.5],
    # 2b
    [1 / 3, 2 / 3, z2_2b],
    [2 / 3, 1 / 3, z2_2b + 0.5],
]
numbers = [0, 0, 1, 1]
cell = moyopy.Cell(basis, positions, numbers)

dataset = moyopy.MoyoDataset(cell, symprec=1e-5, angle_tolerance=None, setting=None)
assert dataset.number == 186
assert dataset.hall_number == 480

hall_symbol_entry = moyopy.HallSymbolEntry(hall_number=dataset.hall_number)
assert hall_symbol_entry.hm_short == "P 6_3 m c"
