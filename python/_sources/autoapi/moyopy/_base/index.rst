moyopy._base
============

.. py:module:: moyopy._base


Classes
-------

.. autoapisummary::

   moyopy._base.Cell
   moyopy._base.Operations


Module Contents
---------------

.. py:class:: Cell(basis: list[list[float]], positions: list[list[float]], numbers: list[int])

   .. py:property:: basis
      :type: list[list[float]]



   .. py:property:: positions
      :type: list[list[float]]



   .. py:property:: numbers
      :type: list[int]



   .. py:method:: serialize_json() -> str


   .. py:method:: deserialize_json(json_str: str) -> Cell
      :classmethod:



.. py:class:: Operations

   .. py:property:: rotations
      :type: list[list[list[float]]]



   .. py:property:: translations
      :type: list[list[float]]



   .. py:property:: num_operations
      :type: int



   .. py:method:: __len__() -> int


