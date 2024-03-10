from __future__ import annotations

from typing import Literal

###############################################################################
# base
###############################################################################

class Structure:
    def __init__(
        self,
        basis: list[list[float]],
        positions: list[list[float]],
        numbers: list[int],
    ): ...
    @property
    def basis(self) -> list[list[float]]: ...
    @property
    def positions(self) -> list[list[float]]: ...
    @property
    def numbers(self) -> list[int]: ...

class Operations:
    @property
    def rotations(self) -> list[list[list[float]]]: ...
    @property
    def translations(self) -> list[list[float]]: ...

###############################################################################
# data
###############################################################################

Setting = int | Literal["spglib", "standard"]

###############################################################################
# lib
###############################################################################

class MoyoDataset:
    def __init__(
        self,
        structure: Structure,
        symprec: float = 1e-5,
        angle_tolerance: float | None = None,
        setting: Setting | None = None,
    ): ...
    @property
    def number(self) -> int: ...
    @property
    def hall_number(self) -> int: ...
    @property
    def operations(self) -> Operations: ...
    @property
    def orbits(self) -> list[int]: ...
    @property
    def wyckoffs(self) -> list[str]: ...
    @property
    def site_symmetry_symbols(self) -> list[str]: ...
    @property
    def std_cell(self) -> Structure: ...
    @property
    def std_linear(self) -> list[list[float]]: ...
    @property
    def std_origin_shift(self) -> list[float]: ...
    @property
    def std_rotation_matrix(self) -> list[list[float]]: ...
    @property
    def prim_std_cell(self) -> Structure: ...
    @property
    def prim_std_linear(self) -> list[list[float]]: ...
    @property
    def prim_std_origin_shift(self) -> list[float]: ...
    @property
    def mapping_std_prim(self) -> list[int]: ...
