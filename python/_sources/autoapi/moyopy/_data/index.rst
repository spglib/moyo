moyopy._data
============

.. py:module:: moyopy._data


Classes
-------

.. autoapisummary::

   moyopy._data.Setting
   moyopy._data.Centering
   moyopy._data.HallSymbolEntry
   moyopy._data.SpaceGroupType


Functions
---------

.. autoapisummary::

   moyopy._data.operations_from_number


Module Contents
---------------

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


.. py:function:: operations_from_number(number: int, setting: Setting) -> moyopy._base.Operations

