from typing import Any

class Cell:
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
    @property
    def num_atoms(self) -> int: ...
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    @classmethod
    def deserialize_json(cls, json_str: str) -> Cell:
        """Deserialize a JSON string to an object"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Cell:
        """Create an object from a dictionary"""

class CollinearMagneticCell:
    def __init__(
        self,
        basis: list[list[float]],
        positions: list[list[float]],
        numbers: list[int],
        magnetic_moments: list[float],
    ): ...
    @property
    def basis(self) -> list[list[float]]: ...
    @property
    def positions(self) -> list[list[float]]: ...
    @property
    def numbers(self) -> list[int]: ...
    @property
    def magnetic_moments(self) -> list[float]: ...
    @property
    def num_atoms(self) -> int: ...
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    @classmethod
    def deserialize_json(cls, json_str: str) -> CollinearMagneticCell:
        """Deserialize a JSON string to an object"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> CollinearMagneticCell:
        """Create an object from a dictionary"""

class NonCollinearMagneticCell:
    def __init__(
        self,
        basis: list[list[float]],
        positions: list[list[float]],
        numbers: list[int],
        magnetic_moments: list[list[float]],
    ): ...
    @property
    def basis(self) -> list[list[float]]: ...
    @property
    def positions(self) -> list[list[float]]: ...
    @property
    def numbers(self) -> list[int]: ...
    @property
    def magnetic_moments(self) -> list[list[float]]: ...
    @property
    def num_atoms(self) -> int: ...
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    @classmethod
    def deserialize_json(cls, json_str: str) -> NonCollinearMagneticCell:
        """Deserialize a JSON string to an object"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> NonCollinearMagneticCell:
        """Create an object from a dictionary"""

class Operations:
    @property
    def rotations(self) -> list[list[list[float]]]: ...
    @property
    def translations(self) -> list[list[float]]: ...
    @property
    def num_operations(self) -> int: ...
    def __len__(self) -> int: ...
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    @classmethod
    def deserialize_json(cls, json_str: str) -> Operations:
        """Deserialize a JSON string to an object"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Operations:
        """Create an object from a dictionary"""

class MagneticOperations:
    @property
    def rotations(self) -> list[list[list[float]]]: ...
    @property
    def translations(self) -> list[list[float]]: ...
    @property
    def time_reversals(self) -> list[bool]: ...
    @property
    def num_operations(self) -> int: ...
    def __len__(self) -> int: ...
    # Serialization and deserialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    @classmethod
    def deserialize_json(cls, json_str: str) -> MagneticOperations:
        """Deserialize a JSON string to an object"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> MagneticOperations:
        """Create an object from a dictionary"""
