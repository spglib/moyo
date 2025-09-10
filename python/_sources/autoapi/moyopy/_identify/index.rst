moyopy._identify
================

.. py:module:: moyopy._identify


Classes
-------

.. autoapisummary::

   moyopy._identify.PointGroup
   moyopy._identify.SpaceGroup


Module Contents
---------------

.. py:class:: PointGroup(prim_rotations: list[list[int]], *, basis: list[list[float]] | None = None)

   .. py:property:: arithmetic_number
      :type: int



   .. py:property:: prim_trans_mat
      :type: list[list[int]]



.. py:class:: SpaceGroup(prim_rotations: list[list[int]], prim_translations: list[list[float]], *, basis: list[list[float]] | None = None, setting: moyopy._data.Setting | None = None, epsilon: float = 0.0001)

   .. py:property:: number
      :type: int



   .. py:property:: hall_number
      :type: int



   .. py:property:: linear
      :type: list[list[int]]



   .. py:property:: origin_shift
      :type: list[float]



