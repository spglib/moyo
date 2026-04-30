from typing import Any

class Cell:
    """A crystal structure: a lattice plus fractional positions and atomic numbers."""
    def __init__(
        self,
        basis: list[list[float]],
        positions: list[list[float]],
        numbers: list[int],
    ):
        """Create a new ``Cell``.

        Parameters
        ----------
        basis : list[list[float]]
            Row-wise basis vectors of the lattice. ``basis[i]`` is the i-th basis vector.
        positions : list[list[float]]
            Fractional coordinates of each site.
        numbers : list[int]
            Atomic number of each site. Must have the same length as ``positions``.
        """
    @property
    def basis(self) -> list[list[float]]:
        """Row-wise basis vectors of the lattice."""
    @property
    def positions(self) -> list[list[float]]:
        """Fractional coordinates of each site."""
    @property
    def numbers(self) -> list[int]:
        """Atomic number of each site."""
    @property
    def num_atoms(self) -> int:
        """Number of atoms in the cell."""
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    @classmethod
    def deserialize_json(cls, json_str: str) -> Cell:
        """Deserialize an object from a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Cell:
        """Create an object from a dictionary."""

class CollinearMagneticCell:
    """A crystal structure with collinear magnetic moments."""
    def __init__(
        self,
        basis: list[list[float]],
        positions: list[list[float]],
        numbers: list[int],
        magnetic_moments: list[float],
    ):
        """Create a new ``CollinearMagneticCell``.

        Parameters
        ----------
        basis : list[list[float]]
            Row-wise basis vectors of the lattice.
        positions : list[list[float]]
            Fractional coordinates of each site.
        numbers : list[int]
            Atomic number of each site.
        magnetic_moments : list[float]
            Scalar magnetic moment of each site (collinear).
        """
    @property
    def basis(self) -> list[list[float]]:
        """Row-wise basis vectors of the lattice."""
    @property
    def positions(self) -> list[list[float]]:
        """Fractional coordinates of each site."""
    @property
    def numbers(self) -> list[int]:
        """Atomic number of each site."""
    @property
    def magnetic_moments(self) -> list[float]:
        """Scalar magnetic moment of each site."""
    @property
    def num_atoms(self) -> int:
        """Number of atoms in the magnetic cell."""
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    @classmethod
    def deserialize_json(cls, json_str: str) -> CollinearMagneticCell:
        """Deserialize an object from a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> CollinearMagneticCell:
        """Create an object from a dictionary."""

class NonCollinearMagneticCell:
    """A crystal structure with non-collinear (vector) magnetic moments."""
    def __init__(
        self,
        basis: list[list[float]],
        positions: list[list[float]],
        numbers: list[int],
        magnetic_moments: list[list[float]],
    ):
        """Create a new ``NonCollinearMagneticCell``.

        Parameters
        ----------
        basis : list[list[float]]
            Row-wise basis vectors of the lattice.
        positions : list[list[float]]
            Fractional coordinates of each site.
        numbers : list[int]
            Atomic number of each site.
        magnetic_moments : list[list[float]]
            Three-component magnetic moment vector of each site.
        """
    @property
    def basis(self) -> list[list[float]]:
        """Row-wise basis vectors of the lattice."""
    @property
    def positions(self) -> list[list[float]]:
        """Fractional coordinates of each site."""
    @property
    def numbers(self) -> list[int]:
        """Atomic number of each site."""
    @property
    def magnetic_moments(self) -> list[list[float]]:
        """Three-component magnetic moment vector of each site."""
    @property
    def num_atoms(self) -> int:
        """Number of atoms in the magnetic cell."""
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    @classmethod
    def deserialize_json(cls, json_str: str) -> NonCollinearMagneticCell:
        """Deserialize an object from a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> NonCollinearMagneticCell:
        """Create an object from a dictionary."""

class Operations:
    """A list of crystallographic symmetry operations (rotation + translation)."""
    @property
    def rotations(self) -> list[list[list[float]]]:
        """Rotation parts of the symmetry operations."""
    @property
    def translations(self) -> list[list[float]]:
        """Translation parts of the symmetry operations (fractional coordinates)."""
    @property
    def num_operations(self) -> int:
        """Number of symmetry operations."""
    def __len__(self) -> int: ...
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    @classmethod
    def deserialize_json(cls, json_str: str) -> Operations:
        """Deserialize an object from a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Operations:
        """Create an object from a dictionary."""

class MagneticOperations:
    """A list of magnetic symmetry operations (rotation + translation + time reversal)."""
    @property
    def rotations(self) -> list[list[list[float]]]:
        """Rotation parts of the magnetic symmetry operations."""
    @property
    def translations(self) -> list[list[float]]:
        """Translation parts of the magnetic symmetry operations (fractional coordinates)."""
    @property
    def time_reversals(self) -> list[bool]:
        """Time-reversal flag for each magnetic symmetry operation."""
    @property
    def num_operations(self) -> int:
        """Number of magnetic symmetry operations."""
    def __len__(self) -> int: ...
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    @classmethod
    def deserialize_json(cls, json_str: str) -> MagneticOperations:
        """Deserialize an object from a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> MagneticOperations:
        """Create an object from a dictionary."""

class UnimodularTransformation:
    """A unimodular transformation: an integer linear part with determinant +-1 plus an
    origin shift."""
    @property
    def linear(self) -> list[list[int]]:
        """Linear part (integer 3x3 matrix with determinant +-1)."""
    @property
    def origin_shift(self) -> list[float]:
        """Origin shift (fractional coordinates)."""
    def serialize_json(self) -> str:
        """Serialize this object to a JSON string."""
    def as_dict(self) -> dict[str, Any]:
        """Convert this object to a dictionary."""
