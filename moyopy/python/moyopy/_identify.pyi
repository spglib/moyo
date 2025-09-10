from moyopy._data import Setting

class PointGroup:
    def __init__(
        self, prim_rotations: list[list[int]], *, basis: list[list[float]] | None = None
    ): ...
    @property
    def arithmetic_number(self) -> int: ...
    @property
    def prim_trans_mat(self) -> list[list[int]]: ...

class SpaceGroup:
    def __init__(
        self,
        prim_rotations: list[list[int]],
        prim_translations: list[list[float]],
        *,
        basis: list[list[float]] | None = None,
        setting: Setting | None = None,
        epsilon: float = 1e-4,
    ): ...
    @property
    def number(self) -> int: ...
    @property
    def hall_number(self) -> int: ...
    @property
    def linear(self) -> list[list[int]]: ...
    @property
    def origin_shift(self) -> list[float]: ...
