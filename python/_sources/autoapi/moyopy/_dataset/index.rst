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

.. py:class:: MoyoDataset(cell: moyopy._base.Cell, *, symprec: float = 0.0001, angle_tolerance: float | None = None, setting: moyopy._data.Setting | None = None, rotate_basis: bool = True)

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

      The input cell is related to the standardized cell by
      ``(std_linear, std_origin_shift)`` and ``std_rotation_matrix``:

      Lattice::

          std_cell.basis.T = std_rotation_matrix @ cell.basis.T @ std_linear

      Fractional positions::

          x_std = np.linalg.inv(std_linear) @ (x_input - std_origin_shift)

      ``std_rotation_matrix`` is a rigid rotation (orthogonal matrix) applied
      only to the Cartesian lattice basis. It does not affect fractional
      coordinates.


   .. py:property:: std_linear
      :type: list[list[float]]


      Linear part of transformation from the input cell to the standardized cell.


   .. py:property:: std_origin_shift
      :type: list[float]


      Origin shift of transformation from the input cell to the standardized cell.


   .. py:property:: std_rotation_matrix
      :type: list[list[float]]


      Rigid rotation (orthogonal matrix) applied to the lattice basis.


   .. py:property:: pearson_symbol
      :type: str


      Pearson symbol for standardized cell.


   .. py:property:: prim_std_cell
      :type: moyopy._base.Cell


      Primitive standardized cell.

      Same transformation convention as the standardized cell above::

          prim_std_cell.basis.T = std_rotation_matrix @ cell.basis.T @ prim_std_linear
          x_prim_std = np.linalg.inv(prim_std_linear) @ (x_input - prim_std_origin_shift)


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

      Serialize an object to a JSON string.



   .. py:method:: deserialize_json(json_str: str) -> Self
      :classmethod:


      Deserialize an object from a JSON string.



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary.



   .. py:method:: from_dict(obj: dict[str, Any]) -> Self
      :classmethod:


      Create an object from a dictionary.



.. py:class:: MoyoCollinearMagneticDataset(magnetic_cell: moyopy._base.CollinearMagneticCell, *, symprec: float = 0.0001, angle_tolerance: float | None = None, mag_symprec: float | None = None, is_axial: bool = False, rotate_basis: bool = True)

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

      The input magnetic cell is related to the standardized magnetic cell by
      ``(std_linear, std_origin_shift)`` and ``std_rotation_matrix``:

      Lattice::

          std_mag_cell.cell.basis.T = std_rotation_matrix @ mag_cell.cell.basis.T @ std_linear

      Fractional positions::

          x_std = np.linalg.inv(std_linear) @ (x_input - std_origin_shift)

      ``std_rotation_matrix`` is a rigid rotation (orthogonal matrix) applied
      only to the Cartesian lattice basis. It does not affect fractional
      coordinates.


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


      Rigid rotation (orthogonal matrix) applied to the lattice basis.


   .. py:property:: prim_std_mag_cell
      :type: moyopy._base.CollinearMagneticCell


      Primitive standardized magnetic cell.

      Same transformation convention as the standardized magnetic cell above::

          prim_std_mag_cell.cell.basis.T = (
              std_rotation_matrix @ mag_cell.cell.basis.T @ prim_std_linear
          )
          x_prim_std = np.linalg.inv(prim_std_linear) @ (x_input - prim_std_origin_shift)


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

      Serialize an object to a JSON string.



   .. py:method:: deserialize_json(json_str: str) -> Self
      :classmethod:


      Deserialize an object from a JSON string.



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary.



   .. py:method:: from_dict(obj: dict[str, Any]) -> Self
      :classmethod:


      Create an object from a dictionary.



.. py:class:: MoyoNonCollinearMagneticDataset(magnetic_cell: moyopy._base.NonCollinearMagneticCell, *, symprec: float = 0.0001, angle_tolerance: float | None = None, mag_symprec: float | None = None, is_axial: bool = True, rotate_basis: bool = True)

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

      The input magnetic cell is related to the standardized magnetic cell by
      ``(std_linear, std_origin_shift)`` and ``std_rotation_matrix``:

      Lattice::

          std_mag_cell.cell.basis.T = std_rotation_matrix @ mag_cell.cell.basis.T @ std_linear

      Fractional positions::

          x_std = np.linalg.inv(std_linear) @ (x_input - std_origin_shift)

      ``std_rotation_matrix`` is a rigid rotation (orthogonal matrix) applied
      only to the Cartesian lattice basis. It does not affect fractional
      coordinates.


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


      Rigid rotation (orthogonal matrix) applied to the lattice basis.


   .. py:property:: prim_std_mag_cell
      :type: moyopy._base.NonCollinearMagneticCell


      Primitive standardized magnetic cell.

      Same transformation convention as the standardized magnetic cell above::

          prim_std_mag_cell.cell.basis.T = (
              std_rotation_matrix @ mag_cell.cell.basis.T @ prim_std_linear
          )
          x_prim_std = np.linalg.inv(prim_std_linear) @ (x_input - prim_std_origin_shift)


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

      Serialize an object to a JSON string.



   .. py:method:: deserialize_json(json_str: str) -> Self
      :classmethod:


      Deserialize an object from a JSON string.



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary.



   .. py:method:: from_dict(obj: dict[str, Any]) -> Self
      :classmethod:


      Create an object from a dictionary.



