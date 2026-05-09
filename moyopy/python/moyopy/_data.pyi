from typing import Any

from moyopy._base import MagneticOperations, Operations

###############################################################################
# Hall symbol data
###############################################################################

class Setting:
    """Preference for the setting of the space group."""
    @classmethod
    def spglib(cls) -> Setting:
        """The setting of the smallest Hall number."""
    @classmethod
    def standard(cls) -> Setting:
        """Unique axis b, cell choice 1 for monoclinic, hexagonal axes for rhombohedral,
        and origin choice 2 for centrosymmetric space groups."""
    @classmethod
    def hall_number(cls, hall_number: int) -> Setting:
        """Specific Hall number from 1 to 530."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class LayerSetting:
    """Preference for the Hall setting of a layer group."""
    @classmethod
    def spglib(cls) -> LayerSetting:
        """The setting with the smallest layer Hall number for each layer group."""
    @classmethod
    def standard(cls) -> LayerSetting:
        """BCS / ITE standard setting per de la Flor et al., Acta Cryst. A77, 559-571 (2021)."""
    @classmethod
    def hall_number(cls, hall_number: int) -> LayerSetting:
        """Specific layer Hall number from 1 to 116."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class Centering:
    @property
    def order(self) -> int:
        """Order of the centering."""
    @property
    def linear(self) -> list[list[int]]:
        """Transformation matrix."""
    @property
    def lattice_points(self) -> list[list[float]]:
        """Unique lattice points."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class LayerCentering:
    """Centering of a layer-group conventional cell. Only ``P`` (primitive) and
    ``C`` (rectangular-centered) occur for layer groups."""
    @property
    def order(self) -> int:
        """Order of the centering."""
    @property
    def linear(self) -> list[list[int]]:
        """Transformation matrix from the primitive cell to the conventional cell.
        The aperiodic axis ``c`` is left untouched."""
    @property
    def lattice_points(self) -> list[list[float]]:
        """Unique lattice points (in fractional coordinates) of the conventional cell.
        The third (``c``) component is always zero because layer-group centerings
        are purely in-plane."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class HallSymbolEntry:
    """An entry containing space-group information for a specified hall_number."""
    def __init__(self, hall_number: int): ...
    @property
    def hall_number(self) -> int:
        """Number for Hall symbols (1 - 530)."""
    @property
    def number(self) -> int:
        """ITA number for space group types (1 - 230)."""
    @property
    def arithmetic_number(self) -> int:
        """Number for arithmetic crystal classes (1 - 73)."""
    @property
    def setting(self) -> Setting:
        """Setting."""
    @property
    def hall_symbol(self) -> str:
        """Hall symbol."""
    @property
    def hm_short(self) -> str:
        """Hermann-Mauguin symbol in short notation."""
    @property
    def hm_full(self) -> str:
        """Hermann-Mauguin symbol in full notation."""
    @property
    def centering(self) -> Centering:
        """Centering."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class LayerHallSymbolEntry:
    """An entry containing layer-group information for a specified layer ``hall_number``."""
    def __init__(self, hall_number: int): ...
    @property
    def hall_number(self) -> int:
        """Sequential number for layer-group Hall settings (1 - 116)."""
    @property
    def number(self) -> int:
        """Layer-group number (1 - 80)."""
    @property
    def arithmetic_number(self) -> int:
        """Number for layer arithmetic crystal classes (1 - 43)."""
    @property
    def setting(self) -> str:
        """Setting code (paper Table 5 axis/origin labels: ``""``, ``"a"``, ``"b"``,
        ``"b-ac"``, ``"c"``, ``"c1"``, ``"c2"``, ``"c3"``, ``"1"``, ``"2"``)."""
    @property
    def hall_symbol(self) -> str:
        """Layer Hall symbol with lowercase ``p``/``c`` lattice prefix."""
    @property
    def hm_short(self) -> str:
        """Hermann-Mauguin symbol in short notation."""
    @property
    def hm_full(self) -> str:
        """Hermann-Mauguin symbol in full notation."""
    @property
    def centering(self) -> LayerCentering:
        """Layer centering."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

###############################################################################
# Group data
###############################################################################

class SpaceGroupType:
    """Space-group type information."""
    def __init__(self, number: int): ...
    # Space-group type
    @property
    def number(self) -> int:
        """ITA number for space group types (1 - 230)."""
    @property
    def hm_short(self) -> str:
        """Hermann-Mauguin symbol in short notation."""
    @property
    def hm_full(self) -> str:
        """Hermann-Mauguin symbol in full notation."""
    # Arithmetic crystal class
    @property
    def arithmetic_number(self) -> int:
        """Number for arithmetic crystal classes (1 - 73)."""
    @property
    def arithmetic_symbol(self) -> str:
        """Symbol for arithmetic crystal class.

        See https://github.com/spglib/moyo/blob/main/moyo/src/data/arithmetic_crystal_class.rs
        for string values.
        """
    # Other classifications
    @property
    def geometric_crystal_class(self) -> str:
        """Geometric crystal class.

        See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
        for string values.
        """
    @property
    def crystal_system(self) -> str:
        """Crystal system.

        See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
        for string values.
        """
    @property
    def bravais_class(self) -> str:
        """Bravais class.

        See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
        for string values.
        """
    @property
    def lattice_system(self) -> str:
        """Lattice system.

        See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
        for string values.
        """
    @property
    def crystal_family(self) -> str:
        """Crystal family.

        See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
        for string values.
        """
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class LayerGroupType:
    """Layer-group type information."""
    def __init__(self, number: int): ...
    # Layer-group type
    @property
    def number(self) -> int:
        """Layer-group number (1 - 80)."""
    @property
    def hm_short(self) -> str:
        """Hermann-Mauguin symbol in short notation."""
    @property
    def hm_full(self) -> str:
        """Hermann-Mauguin symbol in full notation."""
    # Layer arithmetic crystal class
    @property
    def arithmetic_number(self) -> int:
        """Number for layer arithmetic crystal classes (1 - 43)."""
    @property
    def arithmetic_symbol(self) -> str:
        """Symbol for the layer arithmetic crystal class."""
    # Other classifications
    @property
    def geometric_crystal_class(self) -> str:
        """Geometric crystal class. Cubic classes never occur for layer groups.

        See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
        for string values.
        """
    @property
    def layer_bravais_class(self) -> str:
        """Bravais class for the layer group's 2D lattice (one of ``"mp"``,
        ``"op"``, ``"oc"``, ``"tp"``, ``"hp"``)."""
    @property
    def layer_lattice_system(self) -> str:
        """Layer lattice system (one of ``"Oblique"``, ``"Rectangular"``,
        ``"Square"``, ``"Hexagonal"``)."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class MagneticSpaceGroupType:
    """Magnetic space-group type information."""
    def __init__(self, uni_number: int): ...
    @property
    def uni_number(self) -> int:
        """Serial number of UNI (and BNS) symbols."""
    @property
    def litvin_number(self) -> int:
        """Serial number in Litvin's `Magnetic group tables <https://www.iucr.org/publ/978-0-9553602-2-0>`_."""
    @property
    def bns_number(self) -> str:
        """BNS number e.g. '151.32'"""
    @property
    def og_number(self) -> str:
        """OG number e.g. '153.4.1270'"""
    @property
    def number(self) -> int:
        """ITA number for reference space group in BNS setting."""
    @property
    def construct_type(self) -> int:
        """Construct type of magnetic space group from 1 to 4."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class ArithmeticCrystalClass:
    """Arithmetic crystal class information."""
    def __init__(self, arithmetic_number: int): ...
    @property
    def arithmetic_number(self) -> int:
        """Number for arithmetic crystal classes (1 - 73)."""
    @property
    def arithmetic_symbol(self) -> str:
        """Symbol for arithmetic crystal class."""
    @property
    def geometric_crystal_class(self) -> str:
        """Geometric crystal class."""
    @property
    def bravais_class(self) -> str:
        """Bravais class."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

class LayerArithmeticCrystalClass:
    """Layer arithmetic crystal class information."""
    def __init__(self, arithmetic_number: int): ...
    @property
    def arithmetic_number(self) -> int:
        """Number for layer arithmetic crystal classes (1 - 43)."""
    @property
    def symbol(self) -> str:
        """Symbol for the layer arithmetic crystal class (e.g. ``"p1"``,
        ``"c2/m11"``)."""
    @property
    def geometric_crystal_class(self) -> str:
        """Geometric crystal class. Cubic classes never occur for layer groups."""
    @property
    def layer_bravais_class(self) -> str:
        """Bravais class for the layer group's 2D lattice (one of ``"mp"``,
        ``"op"``, ``"oc"``, ``"tp"``, ``"hp"``)."""
    @property
    def layer_lattice_system(self) -> str:
        """Layer lattice system (one of ``"Oblique"``, ``"Rectangular"``,
        ``"Square"``, ``"Hexagonal"``)."""
    # Serialization
    def serialize_json(self) -> str:
        """Serialize an object to a JSON string"""
    def as_dict(self) -> dict[str, Any]:
        """Convert an object to a dictionary"""

###############################################################################
# Misc
###############################################################################

def operations_from_number(
    number: int, *, setting: Setting | None = None, primitive: bool = False
) -> Operations: ...
def operations_from_layer_number(
    number: int, *, setting: LayerSetting | None = None, primitive: bool = False
) -> Operations: ...
def magnetic_operations_from_uni_number(
    uni_number: int, *, primitive: bool = False
) -> MagneticOperations: ...
