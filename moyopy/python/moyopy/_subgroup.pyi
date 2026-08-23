from moyopy._base import Operations

class KlassengleicheSubgroup:
    """A Klassengleiche subgroup embedded in its parent space group."""
    @property
    def operations(self) -> Operations:
        """Coset representatives in the subgroup primitive basis."""
    @property
    def parent_operations(self) -> Operations:
        """The same symmetry operations expressed in the parent primitive basis."""
    @property
    def operation_indices(self) -> list[int]:
        """Mapping from ``operations`` to the input parent operations."""
    @property
    def transformation(self) -> list[list[int]]:
        """Basis whose columns generate the translation sublattice."""
    @property
    def klassengleiche_index(self) -> int:
        """Klassengleiche index, equal to the translation-sublattice index."""

class KlassengleicheSubgroupConjugate:
    """A conjugate embedded Klassengleiche subgroup."""
    @property
    def subgroup(self) -> KlassengleicheSubgroup:
        """The conjugate subgroup embedded in the input parent."""
    @property
    def conjugator(self) -> Operations:
        """Single parent operation whose inverse gives this conjugate."""

class KlassengleicheSubgroupConjugacyClass:
    """A parent-conjugacy class of embedded Klassengleiche subgroups."""
    @property
    def representative(self) -> KlassengleicheSubgroup:
        """Deterministically selected representative of this conjugacy class."""
    @property
    def conjugates(self) -> list[KlassengleicheSubgroupConjugate]:
        """Every distinct conjugate embedding, with the representative first."""
    @property
    def normalizer_operations(self) -> Operations:
        """Normalizer representatives modulo the subgroup translation lattice."""
    @property
    def normalizer_permutations(self) -> list[list[int]]:
        """Induced permutations aligned with ``normalizer_operations``."""

def enumerate_klassengleiche_subgroups(
    prim_rotations: list[list[list[int]]],
    prim_translations: list[list[float]],
    transformation: list[list[int]],
    *,
    epsilon: float = 1e-4,
) -> list[KlassengleicheSubgroupConjugacyClass]:
    """Enumerate Klassengleiche subgroups for a translation sublattice."""

def enumerate_klassengleiche_subgroups_by_index(
    prim_rotations: list[list[list[int]]],
    prim_translations: list[list[float]],
    index: int,
    *,
    epsilon: float = 1e-4,
) -> list[KlassengleicheSubgroupConjugacyClass]:
    """Enumerate Klassengleiche subgroups for every sublattice of an index."""

class TranslationengleicheSubgroup:
    """A Translationengleiche subgroup embedded in its parent space group."""
    @property
    def operations(self) -> Operations:
        """Coset representatives of the subgroup, in the parent primitive basis."""
    @property
    def operation_indices(self) -> list[int]:
        """Mapping from ``operations`` to the input parent operations."""
    @property
    def translationengleiche_index(self) -> int:
        """Translationengleiche index ``[G/T : H/T]``."""
    @property
    def depth(self) -> int:
        """Shortest covering-chain length from the parent in the subgroup lattice."""

class TranslationengleicheSubgroupConjugate:
    """A conjugate embedded Translationengleiche subgroup."""
    @property
    def subgroup(self) -> TranslationengleicheSubgroup:
        """The conjugate subgroup embedded in the input parent."""
    @property
    def conjugator_index(self) -> int:
        """Parent operation whose inverse conjugates the representative here."""

class TranslationengleicheSubgroupConjugacyClass:
    """A parent-conjugacy class of embedded Translationengleiche subgroups."""
    @property
    def representative(self) -> TranslationengleicheSubgroup:
        """Deterministically selected representative of this conjugacy class."""
    @property
    def conjugates(self) -> list[TranslationengleicheSubgroupConjugate]:
        """Every distinct conjugate embedding, with the representative first."""
    @property
    def normalizer_indices(self) -> list[int]:
        """Parent operations that normalize the representative."""
    @property
    def normalizer_permutations(self) -> list[list[int]]:
        """Induced permutations aligned with ``normalizer_indices``."""

def enumerate_translationengleiche_subgroups(
    prim_rotations: list[list[list[int]]],
    prim_translations: list[list[float]],
    *,
    epsilon: float = 1e-4,
) -> list[TranslationengleicheSubgroupConjugacyClass]:
    """Enumerate embedded Translationengleiche subgroups by parent conjugacy."""
