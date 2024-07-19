from __future__ import annotations

import numpy as np
from pymatgen.core import Structure

import moyopy
from moyopy.interface import MoyoAdapter


class MoyoSpacegroupAnalyzer:
    def __init__(
        self,
        structure: Structure,
        symprec: float = 1e-5,
        angle_tolerance: float | None = None,
        setting: moyopy.Setting | None = None,
    ):
        self._cell = MoyoAdapter.from_structure(structure)
        self._dataset = moyopy.MoyoDataset(
            cell=self._cell,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            setting=setting,
        )

    @property
    def number(self) -> int:
        return self._dataset.number

    @property
    def hall_number(self) -> int:
        return self._dataset.hall_number

    @property
    def rotations(self) -> np.ndarray:
        return np.array(self._dataset.operations.rotations)

    @property
    def translations(self) -> np.ndarray:
        return np.array(self._dataset.operations.translations)

    @property
    def orbits(self) -> list[int]:
        return self._dataset.orbits

    @property
    def wyckoffs(self) -> list[str]:
        return self._dataset.wyckoffs

    @property
    def site_symmetry_symbols(self) -> list[str]:
        return self._dataset.site_symmetry_symbols

    @property
    def std_structure(self) -> Structure:
        return MoyoAdapter.get_structure(self._dataset.std_cell)

    @property
    def std_linear(self) -> np.ndarray:
        return np.array(self._dataset.std_linear)

    @property
    def std_origin_shift(self) -> np.ndarray:
        return np.array(self._dataset.std_origin_shift)

    @property
    def std_rotation_matrix(self) -> np.ndarray:
        return np.array(self._dataset.std_rotation_matrix)

    @property
    def prim_std_structure(self) -> Structure:
        return MoyoAdapter.get_structure(self._dataset.prim_std_cell)

    @property
    def prim_std_linear(self) -> np.ndarray:
        return np.array(self._dataset.prim_std_linear)

    @property
    def prim_std_origin_shift(self) -> np.ndarray:
        return np.array(self._dataset.prim_std_origin_shift)

    @property
    def mapping_std_prim(self) -> list[int]:
        return self._dataset.mapping_std_prim

    @property
    def symprec(self) -> float:
        return self._dataset.symprec

    @property
    def angle_tolerance(self) -> float | None:
        return self._dataset.angle_tolerance


if __name__ == "__main__":
    a = 4.0
    structure = Structure(
        lattice=np.eye(3) * a,
        species=["Al"] * 4,
        coords=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ]
        ),
    )
    msa = MoyoSpacegroupAnalyzer(structure)
    assert msa.number == 225
    assert msa.rotations.shape == (48 * 4, 3, 3)
    assert msa.orbits == [0, 0, 0, 0]
    assert msa.prim_std_structure.num_sites == 1
