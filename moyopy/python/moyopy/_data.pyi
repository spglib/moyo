from moyopy._base import Operations

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

###############################################################################
# Group data
###############################################################################

class SpaceGroupType:
    """Space-group type information."""
    def __init__(self, number: int): ...
    # Space group type
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

###############################################################################
# Misc
###############################################################################

def operations_from_number(
    number: int, *, setting: Setting | None = None, primitive: bool = False
) -> Operations: ...
