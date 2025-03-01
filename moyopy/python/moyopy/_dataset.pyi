from typing import Any

from typing_extensions import Self

from moyopy._base import (
    Cell,
    CollinearMagneticCell,
    MagneticOperations,
    NonCollinearMagneticCell,
    Operations,
)
from moyopy._data import (
    Setting,
)

class MoyoDataset:
    """A dataset containing symmetry information of the input crystal structure."""
    def __init__(
        self,
        cell: Cell,
        *,
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
            Symmetry search tolerance in the unit of cell.basis.
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
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize the dataset to a JSON string."""
    @classmethod
    def deserialize_json(cls, json_str: str) -> Self:
        """Deserialize the dataset from a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert the dataset to a dictionary."""
    @classmethod
    def from_dict(cls, obj: dict[str, Any]) -> Self:
        """Create a dataset from a dictionary."""

class MoyoCollinearMagneticDataset:
    """A dataset containing magnetic symmetry information of the input collinear magnetic
    structure."""
    def __init__(
        self,
        magnetic_cell: CollinearMagneticCell,
        *,
        symprec: float = 1e-4,
        angle_tolerance: float | None = None,
        mag_symprec: float | None = None,
        is_axial: bool = False,
    ):
        """
        Parameters
        ----------
        magnetic_cell: CollinearMagneticCell
            Input collinear magnetic structure.
        symprec: float
            Symmetry search tolerance in the unit of magnetic_cell.basis.
        angle_tolerance: float | None
            Symmetry search tolerance in the unit of radians.
        mag_symprec: float | None
            Symmetry search tolerance in the unit of magnetic moments.
        is_axial: bool
            Whether the magnetic moments are axial on improper operations.
        """
    # Magnetic space-group type
    @property
    def uni_number(self) -> int:
        """UNI number for magnetic space-group type."""
    # Magnetic symmetry operations in the input cell
    @property
    def magnetic_operations(self) -> MagneticOperations:
        """Magnetic symmetry operations in the input cell."""
    # Site symmetry
    @property
    def orbits(self) -> list[int]:
        """The `i`th atom in the input magnetic cell is equivalent to the `orbits[i]`th atom
        in the **input** magnetic cell. For example, orbits=[0, 0, 2, 2, 2, 2] means
        the first two atoms are equivalent and the last four atoms are equivalent to each other.
        """
    # Standardized magnetic cell
    @property
    def std_mag_cell(self) -> CollinearMagneticCell:
        """Standardized magnetic cell."""
    @property
    def std_linear(self) -> list[list[float]]:
        """Linear part of transformation from the input magnetic cell to the standardized
        magnetic cell."""
    @property
    def std_origin_shift(self) -> list[float]:
        """Origin shift of transformation from the input magnetic cell to the standardized
        magnetic cell."""
    @property
    def std_rotation_matrix(self) -> list[list[float]]:
        """Rigid rotation."""
    # Primitive standardized magnetic cell
    @property
    def prim_std_mag_cell(self) -> CollinearMagneticCell:
        """Primitive standardized magnetic cell."""
    @property
    def prim_std_linear(self) -> list[list[float]]:
        """Linear part of transformation from the input magnetic cell to the primitive
        standardized magnetic cell."""
    @property
    def prim_std_origin_shift(self) -> list[float]:
        """Origin shift of transformation from the input magnetic cell to the primitive
        standardized magnetic cell."""
    @property
    def mapping_std_prim(self) -> list[int]:
        """Mapping sites in the input magnetic cell to those in the primitive standardized magnetic
        cell. The `i`th atom in the input magnetic cell is mapped to the `mapping_to_std_prim[i]`th
        atom in the primitive standardized magnetic cell.
        """
    # Final parameters
    @property
    def symprec(self) -> float:
        """Actually used `symprec` in iterative symmetry search."""
    @property
    def angle_tolerance(self) -> float | None:
        """Actually used `angle_tolerance` in iterative symmetry search."""
    @property
    def mag_symprec(self) -> float | None:
        """Actually used `mag_symprec` in iterative symmetry search."""
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize the dataset to a JSON string."""
    @classmethod
    def deserialize_json(cls, json_str: str) -> Self:
        """Deserialize the dataset from a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert the dataset to a dictionary."""
    @classmethod
    def from_dict(cls, obj: dict[str, Any]) -> Self:
        """Create a dataset from a dictionary."""

class MoyoNonCollinearMagneticDataset:
    """A dataset containing magnetic symmetry information of the input non-collinear magnetic
    structure."""
    def __init__(
        self,
        magnetic_cell: NonCollinearMagneticCell,
        *,
        symprec: float = 1e-4,
        angle_tolerance: float | None = None,
        mag_symprec: float | None = None,
        is_axial: bool = True,
    ):
        """
        Parameters
        ----------
        magnetic_cell: NonCollinearMagneticCell
            Input non-collinear magnetic structure.
        symprec: float
            Symmetry search tolerance in the unit of magnetic_cell.basis.
        angle_tolerance: float | None
            Symmetry search tolerance in the unit of radians.
        mag_symprec: float | None
            Symmetry search tolerance in the unit of magnetic moments.
        is_axial: bool
            Whether the magnetic moments are axial on improper operations.
        """
    # Magnetic space-group type
    @property
    def uni_number(self) -> int:
        """UNI number for magnetic space-group type."""
    # Magnetic symmetry operations in the input cell
    @property
    def magnetic_operations(self) -> MagneticOperations:
        """Magnetic symmetry operations in the input cell."""
    # Site symmetry
    @property
    def orbits(self) -> list[int]:
        """The `i`th atom in the input magnetic cell is equivalent to the `orbits[i]`th atom
        in the **input** magnetic cell. For example, orbits=[0, 0, 2, 2, 2, 2] means
        the first two atoms are equivalent and the last four atoms are equivalent to each other.
        """
    # Standardized magnetic cell
    @property
    def std_mag_cell(self) -> NonCollinearMagneticCell:
        """Standardized magnetic cell."""
    @property
    def std_linear(self) -> list[list[float]]:
        """Linear part of transformation from the input magnetic cell to the standardized
        magnetic cell."""
    @property
    def std_origin_shift(self) -> list[float]:
        """Origin shift of transformation from the input magnetic cell to the standardized
        magnetic cell."""
    @property
    def std_rotation_matrix(self) -> list[list[float]]:
        """Rigid rotation."""
    # Primitive standardized magnetic cell
    @property
    def prim_std_mag_cell(self) -> NonCollinearMagneticCell:
        """Primitive standardized magnetic cell."""
    @property
    def prim_std_linear(self) -> list[list[float]]:
        """Linear part of transformation from the input magnetic cell to the primitive
        standardized magnetic cell."""
    @property
    def prim_std_origin_shift(self) -> list[float]:
        """Origin shift of transformation from the input magnetic cell to the primitive
        standardized magnetic cell."""
    @property
    def mapping_std_prim(self) -> list[int]:
        """Mapping sites in the input magnetic cell to those in the primitive standardized magnetic
        cell. The `i`th atom in the input magnetic cell is mapped to the `mapping_to_std_prim[i]`th
        atom in the primitive standardized magnetic cell.
        """
    # Final parameters
    @property
    def symprec(self) -> float:
        """Actually used `symprec` in iterative symmetry search."""
    @property
    def angle_tolerance(self) -> float | None:
        """Actually used `angle_tolerance` in iterative symmetry search."""
    @property
    def mag_symprec(self) -> float | None:
        """Actually used `mag_symprec` in iterative symmetry search."""
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize the dataset to a JSON string."""
    @classmethod
    def deserialize_json(cls, json_str: str) -> Self:
        """Deserialize the dataset from a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert the dataset to a dictionary."""
    @classmethod
    def from_dict(cls, obj: dict[str, Any]) -> Self:
        """Create a dataset from a dictionary."""
