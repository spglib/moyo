moyopy._dataset
===============

.. py:module:: moyopy._dataset


Classes
-------

.. autoapisummary::

   moyopy._dataset.MoyoDataset
   moyopy._dataset.MoyoCollinearMagneticDataset
   moyopy._dataset.MoyoNonCollinearMagneticDataset


Module Contents
---------------

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



