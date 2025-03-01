moyopy._base
============

.. py:module:: moyopy._base


Classes
-------

.. autoapisummary::

   moyopy._base.Cell
   moyopy._base.CollinearMagneticCell
   moyopy._base.NonCollinearMagneticCell
   moyopy._base.Operations
   moyopy._base.MagneticOperations


Module Contents
---------------

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


.. py:class:: MagneticOperations

   .. py:property:: rotations
      :type: list[list[list[float]]]



   .. py:property:: translations
      :type: list[list[float]]



   .. py:property:: time_reversals
      :type: list[bool]



   .. py:property:: num_operations
      :type: int



   .. py:method:: __len__() -> int


