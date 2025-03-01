moyopy
======

.. py:module:: moyopy


Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/moyopy/_base/index
   /autoapi/moyopy/_data/index
   /autoapi/moyopy/_dataset/index
   /autoapi/moyopy/_moyopy/index
   /autoapi/moyopy/interface/index


Attributes
----------

.. autoapisummary::

   moyopy.__version__


Classes
-------

.. autoapisummary::

   moyopy.Cell
   moyopy.CollinearMagneticCell
   moyopy.NonCollinearMagneticCell
   moyopy.Operations
   moyopy.Setting
   moyopy.Centering
   moyopy.HallSymbolEntry
   moyopy.SpaceGroupType
   moyopy.MagneticSpaceGroupType
   moyopy.MoyoDataset
   moyopy.MoyoCollinearMagneticDataset
   moyopy.MoyoNonCollinearMagneticDataset


Functions
---------

.. autoapisummary::

   moyopy.operations_from_number


Package Contents
----------------

.. py:class:: Cell(basis: list[list[float]], positions: list[list[float]], numbers: list[int])

   .. py:property:: basis
      :type: list[list[float]]



   .. py:property:: positions
      :type: list[list[float]]



   .. py:property:: numbers
      :type: list[int]



   .. py:property:: num_atoms
      :type: int



   .. py:method:: serialize_json() -> str

      Serialize the cell to a JSON string



   .. py:method:: deserialize_json(json_str: str) -> Cell
      :classmethod:


      Deserialize a JSON string to a Cell object



   .. py:method:: as_dict() -> dict[str, Any]

      Convert the cell to a dictionary



   .. py:method:: from_dict(data: dict[str, Any]) -> Cell
      :classmethod:


      Create a cell from a dictionary



.. py:class:: CollinearMagneticCell(basis: list[list[float]], positions: list[list[float]], numbers: list[int], magnetic_moments: list[float])

   .. py:property:: basis
      :type: list[list[float]]



   .. py:property:: positions
      :type: list[list[float]]



   .. py:property:: numbers
      :type: list[int]



   .. py:property:: magnetic_moments
      :type: list[float]



   .. py:property:: num_atoms
      :type: int



   .. py:method:: serialize_json() -> str

      Serialize the cell to a JSON string



   .. py:method:: deserialize_json(json_str: str) -> Cell
      :classmethod:


      Deserialize a JSON string to a Cell object



   .. py:method:: as_dict() -> dict[str, Any]

      Convert the cell to a dictionary



   .. py:method:: from_dict(data: dict[str, Any]) -> Cell
      :classmethod:


      Create a cell from a dictionary



.. py:class:: NonCollinearMagneticCell(basis: list[list[float]], positions: list[list[float]], numbers: list[int], magnetic_moments: list[list[float]])

   .. py:property:: basis
      :type: list[list[float]]



   .. py:property:: positions
      :type: list[list[float]]



   .. py:property:: numbers
      :type: list[int]



   .. py:property:: magnetic_moments
      :type: list[list[float]]



   .. py:property:: num_atoms
      :type: int



   .. py:method:: serialize_json() -> str

      Serialize the cell to a JSON string



   .. py:method:: deserialize_json(json_str: str) -> Cell
      :classmethod:


      Deserialize a JSON string to a Cell object



   .. py:method:: as_dict() -> dict[str, Any]

      Convert the cell to a dictionary



   .. py:method:: from_dict(data: dict[str, Any]) -> Cell
      :classmethod:


      Create a cell from a dictionary



.. py:class:: Operations

   .. py:property:: rotations
      :type: list[list[list[float]]]



   .. py:property:: translations
      :type: list[list[float]]



   .. py:property:: num_operations
      :type: int



   .. py:method:: __len__() -> int


.. py:class:: Setting

   Preference for the setting of the space group.


   .. py:method:: spglib() -> Setting
      :classmethod:


      The setting of the smallest Hall number.



   .. py:method:: standard() -> Setting
      :classmethod:


      Unique axis b, cell choice 1 for monoclinic, hexagonal axes for rhombohedral,
      and origin choice 2 for centrosymmetric space groups.



   .. py:method:: hall_number(hall_number: int) -> Setting
      :classmethod:


      Specific Hall number from 1 to 530.



.. py:class:: Centering

.. py:class:: HallSymbolEntry(hall_number: int)

   An entry containing space-group information for a specified hall_number.


   .. py:property:: hall_number
      :type: int


      Number for Hall symbols (1 - 530).


   .. py:property:: number
      :type: int


      ITA number for space group types (1 - 230).


   .. py:property:: arithmetic_number
      :type: int


      Number for arithmetic crystal classes (1 - 73).


   .. py:property:: setting
      :type: Setting


      Setting.


   .. py:property:: hall_symbol
      :type: str


      Hall symbol.


   .. py:property:: hm_short
      :type: str


      Hermann-Mauguin symbol in short notation.


   .. py:property:: hm_full
      :type: str


      Hermann-Mauguin symbol in full notation.


   .. py:property:: centering
      :type: Centering


      Centering.


.. py:class:: SpaceGroupType(number: int)

   Space-group type information.


   .. py:property:: number
      :type: int


      ITA number for space group types (1 - 230).


   .. py:property:: hm_short
      :type: str


      Hermann-Mauguin symbol in short notation.


   .. py:property:: hm_full
      :type: str


      Hermann-Mauguin symbol in full notation.


   .. py:property:: arithmetic_number
      :type: int


      Number for arithmetic crystal classes (1 - 73).


   .. py:property:: arithmetic_symbol
      :type: str


      Symbol for arithmetic crystal class.

      See https://github.com/spglib/moyo/blob/main/moyo/src/data/arithmetic_crystal_class.rs
      for string values.


   .. py:property:: geometric_crystal_class
      :type: str


      Geometric crystal class.

      See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
      for string values.


   .. py:property:: crystal_system
      :type: str


      Crystal system.

      See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
      for string values.


   .. py:property:: bravais_class
      :type: str


      Bravais class.

      See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
      for string values.


   .. py:property:: lattice_system
      :type: str


      Lattice system.

      See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
      for string values.


   .. py:property:: crystal_family
      :type: str


      Crystal family.

      See https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs
      for string values.


.. py:class:: MagneticSpaceGroupType(uni_number: int)

   Magnetic space-group type information.


   .. py:property:: uni_number
      :type: int


      Serial number of UNI (and BNS) symbols.


   .. py:property:: litvin_number
      :type: int


      Serial number in Litvin's `Magnetic group tables <https://www.iucr.org/publ/978-0-9553602-2-0>`_.


   .. py:property:: bns_number
      :type: str


      BNS number e.g. '151.32'


   .. py:property:: og_number
      :type: str


      OG number e.g. '153.4.1270'


   .. py:property:: number
      :type: int


      ITA number for reference space group in BNS setting.


   .. py:property:: construct_type
      :type: int


      Construct type of magnetic space group from 1 to 4.


.. py:function:: operations_from_number(number: int, *, setting: Setting | None = None) -> moyopy._base.Operations

.. py:class:: MoyoDataset(cell: moyopy._base.Cell, *, symprec: float = 0.0001, angle_tolerance: float | None = None, setting: moyopy._data.Setting | None = None)

   A dataset containing symmetry information of the input crystal structure.


   .. py:property:: number
      :type: int


      Space group number.


   .. py:property:: hall_number
      :type: int


      Hall symbol number.


   .. py:property:: operations
      :type: moyopy._base.Operations


      Symmetry operations in the input cell.


   .. py:property:: orbits
      :type: list[int]


      Spglib's `crystallographic_orbits` not `equivalent_atoms`.

      The `i`th atom in the input cell is equivalent to the `orbits[i]`th atom in the **input**
      cell. For example, orbits=[0, 0, 2, 2, 2, 2] means the first two atoms are equivalent
      and the last four atoms are equivalent to each other.


   .. py:property:: wyckoffs
      :type: list[str]


      Wyckoff letters for each site in the input cell.


   .. py:property:: site_symmetry_symbols
      :type: list[str]


      Site symmetry symbols for each site in the input cell.

      The orientation of the site symmetry is w.r.t. the standardized cell.


   .. py:property:: std_cell
      :type: moyopy._base.Cell


      Standardized cell.


   .. py:property:: std_linear
      :type: list[list[float]]


      Linear part of transformation from the input cell to the standardized cell.


   .. py:property:: std_origin_shift
      :type: list[float]


      Origin shift of transformation from the input cell to the standardized cell.


   .. py:property:: std_rotation_matrix
      :type: list[list[float]]


      Rigid rotation.


   .. py:property:: pearson_symbol
      :type: str


      Pearson symbol for standardized cell.


   .. py:property:: prim_std_cell
      :type: moyopy._base.Cell


      Primitive standardized cell.


   .. py:property:: prim_std_linear
      :type: list[list[float]]


      Linear part of transformation from the input cell to the primitive standardized cell.


   .. py:property:: prim_std_origin_shift
      :type: list[float]


      Origin shift of transformation from the input cell to the primitive standardized
      cell.


   .. py:property:: mapping_std_prim
      :type: list[int]


      Mapping sites in the input cell to those in the primitive standardized cell.

      The `i`th atom in the input cell is mapped to the `mapping_to_std_prim[i]`th atom in the
      primitive standardized cell.


   .. py:property:: symprec
      :type: float


      Actually used `symprec` in iterative symmetry search.


   .. py:property:: angle_tolerance
      :type: float | None


      Actually used `angle_tolerance` in iterative symmetry search.


   .. py:method:: serialize_json() -> str

      Serialize the dataset to a JSON string.



   .. py:method:: deserialize_json(json_str: str) -> typing_extensions.Self
      :classmethod:


      Deserialize the dataset from a JSON string.



   .. py:method:: as_dict() -> dict[str, Any]

      Convert the dataset to a dictionary.



   .. py:method:: from_dict(obj: dict[str, Any]) -> typing_extensions.Self
      :classmethod:


      Create a dataset from a dictionary.



.. py:class:: MoyoCollinearMagneticDataset(magnetic_cell: moyopy._base.CollinearMagneticCell, *, symprec: float = 0.0001, angle_tolerance: float | None = None, mag_symprec: float | None = None, is_axial: bool = False)

   A dataset containing magnetic symmetry information of the input collinear magnetic
   structure.


   .. py:property:: uni_number
      :type: int


      UNI number for magnetic space-group type.


   .. py:property:: magnetic_operations
      :type: moyopy._base.MagneticOperations


      Magnetic symmetry operations in the input cell.


   .. py:property:: orbits
      :type: list[int]


      The `i`th atom in the input magnetic cell is equivalent to the `orbits[i]`th atom
      in the **input** magnetic cell. For example, orbits=[0, 0, 2, 2, 2, 2] means
      the first two atoms are equivalent and the last four atoms are equivalent to each other.


   .. py:property:: std_mag_cell
      :type: moyopy._base.CollinearMagneticCell


      Standardized magnetic cell.


   .. py:property:: std_linear
      :type: list[list[float]]


      Linear part of transformation from the input magnetic cell to the standardized
      magnetic cell.


   .. py:property:: std_origin_shift
      :type: list[float]


      Origin shift of transformation from the input magnetic cell to the standardized
      magnetic cell.


   .. py:property:: std_rotation_matrix
      :type: list[list[float]]


      Rigid rotation.


   .. py:property:: prim_std_mag_cell
      :type: moyopy._base.CollinearMagneticCell


      Primitive standardized magnetic cell.


   .. py:property:: prim_std_linear
      :type: list[list[float]]


      Linear part of transformation from the input magnetic cell to the primitive
      standardized magnetic cell.


   .. py:property:: prim_std_origin_shift
      :type: list[float]


      Origin shift of transformation from the input magnetic cell to the primitive
      standardized magnetic cell.


   .. py:property:: mapping_std_prim
      :type: list[int]


      Mapping sites in the input magnetic cell to those in the primitive standardized magnetic
      cell. The `i`th atom in the input magnetic cell is mapped to the `mapping_to_std_prim[i]`th
      atom in the primitive standardized magnetic cell.


   .. py:property:: symprec
      :type: float


      Actually used `symprec` in iterative symmetry search.


   .. py:property:: angle_tolerance
      :type: float | None


      Actually used `angle_tolerance` in iterative symmetry search.


   .. py:property:: mag_symprec
      :type: float | None


      Actually used `mag_symprec` in iterative symmetry search.


   .. py:method:: serialize_json() -> str

      Serialize the dataset to a JSON string.



   .. py:method:: deserialize_json(json_str: str) -> typing_extensions.Self
      :classmethod:


      Deserialize the dataset from a JSON string.



   .. py:method:: as_dict() -> dict[str, Any]

      Convert the dataset to a dictionary.



   .. py:method:: from_dict(obj: dict[str, Any]) -> typing_extensions.Self
      :classmethod:


      Create a dataset from a dictionary.



.. py:class:: MoyoNonCollinearMagneticDataset(magnetic_cell: moyopy._base.NonCollinearMagneticCell, *, symprec: float = 0.0001, angle_tolerance: float | None = None, mag_symprec: float | None = None, is_axial: bool = True)

   A dataset containing magnetic symmetry information of the input non-collinear magnetic
   structure.


   .. py:property:: uni_number
      :type: int


      UNI number for magnetic space-group type.


   .. py:property:: magnetic_operations
      :type: moyopy._base.MagneticOperations


      Magnetic symmetry operations in the input cell.


   .. py:property:: orbits
      :type: list[int]


      The `i`th atom in the input magnetic cell is equivalent to the `orbits[i]`th atom
      in the **input** magnetic cell. For example, orbits=[0, 0, 2, 2, 2, 2] means
      the first two atoms are equivalent and the last four atoms are equivalent to each other.


   .. py:property:: std_mag_cell
      :type: moyopy._base.NonCollinearMagneticCell


      Standardized magnetic cell.


   .. py:property:: std_linear
      :type: list[list[float]]


      Linear part of transformation from the input magnetic cell to the standardized
      magnetic cell.


   .. py:property:: std_origin_shift
      :type: list[float]


      Origin shift of transformation from the input magnetic cell to the standardized
      magnetic cell.


   .. py:property:: std_rotation_matrix
      :type: list[list[float]]


      Rigid rotation.


   .. py:property:: prim_std_mag_cell
      :type: moyopy._base.NonCollinearMagneticCell


      Primitive standardized magnetic cell.


   .. py:property:: prim_std_linear
      :type: list[list[float]]


      Linear part of transformation from the input magnetic cell to the primitive
      standardized magnetic cell.


   .. py:property:: prim_std_origin_shift
      :type: list[float]


      Origin shift of transformation from the input magnetic cell to the primitive
      standardized magnetic cell.


   .. py:property:: mapping_std_prim
      :type: list[int]


      Mapping sites in the input magnetic cell to those in the primitive standardized magnetic
      cell. The `i`th atom in the input magnetic cell is mapped to the `mapping_to_std_prim[i]`th
      atom in the primitive standardized magnetic cell.


   .. py:property:: symprec
      :type: float


      Actually used `symprec` in iterative symmetry search.


   .. py:property:: angle_tolerance
      :type: float | None


      Actually used `angle_tolerance` in iterative symmetry search.


   .. py:property:: mag_symprec
      :type: float | None


      Actually used `mag_symprec` in iterative symmetry search.


   .. py:method:: serialize_json() -> str

      Serialize the dataset to a JSON string.



   .. py:method:: deserialize_json(json_str: str) -> typing_extensions.Self
      :classmethod:


      Deserialize the dataset from a JSON string.



   .. py:method:: as_dict() -> dict[str, Any]

      Convert the dataset to a dictionary.



   .. py:method:: from_dict(obj: dict[str, Any]) -> typing_extensions.Self
      :classmethod:


      Create a dataset from a dictionary.



.. py:data:: __version__
   :type:  str

