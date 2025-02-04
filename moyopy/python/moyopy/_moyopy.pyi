from moyopy._base import Cell, Operations  # noqa: F401
from moyopy._data import (
    Centering,
    HallSymbolEntry,
    Setting,
    SpaceGroupType,
    operations_from_number,
)  # noqa: F401

__version__: str

class MoyoDataset:
    """A dataset containing symmetry information of the input crystal structure."""
    def __init__(
        self,
        cell: Cell,
        symprec: float = 1e-4,
        angle_tolerance: float | None = None,
        setting: Setting | None = None,
    ):
        """
        Parameters
        ----------
        cell: Cell
            Input crystal structure.
        symprec: float
            Symmetry search tolerance in the unit of cell.lattice.
        angle_tolerance: float | None
            Symmetry search tolerance in the unit of radians.
        setting: Setting | None
            Preference for the setting of the space group.
        """
    # Space-group type
    @property
    def number(self) -> int:
        """Space group number."""
    @property
    def hall_number(self) -> int:
        """Hall symbol number."""
    # Symmetry operations in the input cell
    @property
    def operations(self) -> Operations:
        """Symmetry operations in the input cell."""
    # Site symmetry
    @property
    def orbits(self) -> list[int]:
        """Spglib's `crystallographic_orbits` not `equivalent_atoms`.

        The `i`th atom in the input cell is equivalent to the `orbits[i]`th atom in the **input**
        cell. For example, orbits=[0, 0, 2, 2, 2, 2] means the first two atoms are equivalent
        and the last four atoms are equivalent to each other.
        """
    @property
    def wyckoffs(self) -> list[str]:
        """Wyckoff letters for each site in the input cell."""
    @property
    def site_symmetry_symbols(self) -> list[str]:
        """Site symmetry symbols for each site in the input cell.

        The orientation of the site symmetry is w.r.t. the standardized cell.
        """
    # Standardized cell
    @property
    def std_cell(self) -> Cell:
        """Standardized cell."""
    @property
    def std_linear(self) -> list[list[float]]:
        """Linear part of transformation from the input cell to the standardized cell."""
    @property
    def std_origin_shift(self) -> list[float]:
        """Origin shift of transformation from the input cell to the standardized cell."""
    @property
    def std_rotation_matrix(self) -> list[list[float]]:
        """Rigid rotation."""
    @property
    def pearson_symbol(self) -> str:
        """Pearson symbol for standardized cell."""
    # Primitive standardized cell
    @property
    def prim_std_cell(self) -> Cell:
        """Primitive standardized cell."""
    @property
    def prim_std_linear(self) -> list[list[float]]:
        """Linear part of transformation from the input cell to the primitive standardized cell."""
    @property
    def prim_std_origin_shift(self) -> list[float]:
        """Origin shift of transformation from the input cell to the primitive standardized
        cell."""
    @property
    def mapping_std_prim(self) -> list[int]:
        """Mapping sites in the input cell to those in the primitive standardized cell.

        The `i`th atom in the input cell is mapped to the `mapping_to_std_prim[i]`th atom in the
        primitive standardized cell.
        """
    # Final parameters
    @property
    def symprec(self) -> float:
        """Actually used `symprec` in iterative symmetry search."""
    @property
    def angle_tolerance(self) -> float | None:
        """Actually used `angle_tolerance` in iterative symmetry search."""

__all__ = [
    # base
    "Cell",
    "Operations",
    # data
    "Setting",
    "Centering",
    "HallSymbolEntry",
    "SpaceGroupType",
    "operations_from_number",
    # lib
    "__version__",
    "MoyoDataset",
]
