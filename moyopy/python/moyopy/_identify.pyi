from typing import Any

from moyopy._base import UnimodularTransformation
from moyopy._data import Setting

class PointGroup:
    """Point group identified from a list of rotation matrices."""
    def __init__(self, prim_rotations: list[list[int]], *, basis: list[list[float]] | None = None):
        """Identify the point group of the given primitive rotations.

        Parameters
        ----------
        prim_rotations : list[list[list[int]]]
            Rotation matrices in the primitive cell.
        basis : list[list[float]] | None
            Row-wise basis vectors of the primitive lattice. If ``None``, an identity basis
            is assumed.
        """
    @property
    def arithmetic_number(self) -> int:
        """Number for the arithmetic crystal class (1 - 73)."""
    @property
    def prim_trans_mat(self) -> list[list[int]]:
        """Transformation matrix from the input primitive basis to the standardized basis."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""

class SpaceGroup:
    """Space group identified from a list of primitive symmetry operations."""
    def __init__(
        self,
        prim_rotations: list[list[int]],
        prim_translations: list[list[float]],
        *,
        basis: list[list[float]] | None = None,
        setting: Setting | None = None,
        epsilon: float = 1e-4,
    ):
        """Identify the space group from primitive rotations and translations.

        Parameters
        ----------
        prim_rotations : list[list[list[int]]]
            Rotation matrices of the symmetry operations of the primitive cell.
        prim_translations : list[list[float]]
            Translation vectors of the symmetry operations of the primitive cell.
        basis : list[list[float]] | None
            Row-wise basis vectors of the primitive lattice. If ``None``, an identity basis
            is assumed.
        setting : Setting | None
            Preference for the standardized setting of the detected space-group type.
        epsilon : float
            Numerical tolerance for matching translations.
        """
    @property
    def number(self) -> int:
        """ITA number for the identified space group (1 - 230)."""
    @property
    def hall_number(self) -> int:
        """Hall symbol number (1 - 530) for the chosen setting."""
    @property
    def linear(self) -> list[list[int]]:
        """Linear part of the transformation from the input primitive basis to the
        standardized basis."""
    @property
    def origin_shift(self) -> list[float]:
        """Origin shift of the transformation from the input primitive basis to the
        standardized basis."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""

class MagneticSpaceGroup:
    """Magnetic space group identified from a list of primitive magnetic operations."""
    def __init__(
        self,
        prim_rotations: list[list[int]],
        prim_translations: list[list[float]],
        prim_time_reversals: list[bool],
        *,
        basis: list[list[float]] | None = None,
        epsilon: float = 1e-4,
    ):
        """Identify the magnetic space group from primitive magnetic operations.

        Parameters
        ----------
        prim_rotations : list[list[list[int]]]
            Rotation matrices of the magnetic operations of the primitive cell.
        prim_translations : list[list[float]]
            Translation vectors of the magnetic operations of the primitive cell.
        prim_time_reversals : list[bool]
            Time-reversal flag for each magnetic operation of the primitive cell.
        basis : list[list[float]] | None
            Row-wise basis vectors of the primitive lattice. If ``None``, an identity basis
            is assumed.
        epsilon : float
            Numerical tolerance for matching translations.
        """
    @property
    def uni_number(self) -> int:
        """Serial number of UNI (and BNS) symbols."""
    @property
    def linear(self) -> list[list[int]]:
        """Linear part of the transformation from the input primitive basis to the
        standardized basis."""
    @property
    def origin_shift(self) -> list[float]:
        """Origin shift of the transformation from the input primitive basis to the
        standardized basis."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""

def integral_normalizer(
    prim_rotations: list[list[list[int]]],
    prim_translations: list[list[float]],
    *,
    prim_generators: list[int] | None = None,
    epsilon: float = 1e-4,
) -> list[UnimodularTransformation]:
    """Compute the integral normalizer of a space group.

    Returns the unimodular transformations that conjugate the given space group into itself.

    Parameters
    ----------
    prim_rotations : list[list[list[int]]]
        Rotation matrices of the symmetry operations of the primitive cell.
    prim_translations : list[list[float]]
        Translation vectors of the symmetry operations of the primitive cell.
    prim_generators : list[int] | None
        Optional indices of operations to use as generators. If ``None``, a small generating
        set is derived automatically.
    epsilon : float
        Numerical tolerance for matching translations.
    """
