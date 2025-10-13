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
   moyopy._data.MagneticSpaceGroupType
   moyopy._data.ArithmeticCrystalClass


Functions
---------

.. autoapisummary::

   moyopy._data.operations_from_number
   moyopy._data.magnetic_operations_from_uni_number


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



   .. py:method:: serialize_json() -> str

      Serialize an object to a JSON string



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary



.. py:class:: Centering

   .. py:property:: order
      :type: int


      Order of the centering.


   .. py:property:: linear
      :type: list[list[int]]


      Transformation matrix.


   .. py:property:: lattice_points
      :type: list[list[float]]


      Unique lattice points.


   .. py:method:: serialize_json() -> str

      Serialize an object to a JSON string



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary



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


   .. py:method:: serialize_json() -> str

      Serialize an object to a JSON string



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary



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


   .. py:method:: serialize_json() -> str

      Serialize an object to a JSON string



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary



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


   .. py:method:: serialize_json() -> str

      Serialize an object to a JSON string



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary



.. py:class:: ArithmeticCrystalClass(arithmetic_number: int)

   Arithmetic crystal class information.


   .. py:property:: arithmetic_number
      :type: int


      Number for arithmetic crystal classes (1 - 73).


   .. py:property:: arithmetic_symbol
      :type: str


      Symbol for arithmetic crystal class.


   .. py:property:: geometric_crystal_class
      :type: str


      Geometric crystal class.


   .. py:property:: bravais_class
      :type: str


      Bravais class.


   .. py:method:: serialize_json() -> str

      Serialize an object to a JSON string



   .. py:method:: as_dict() -> dict[str, Any]

      Convert an object to a dictionary



.. py:function:: operations_from_number(number: int, *, setting: Setting | None = None, primitive: bool = False) -> moyopy._base.Operations

.. py:function:: magnetic_operations_from_uni_number(uni_number: int, *, primitive: bool = False) -> moyopy._base.MagneticOperations

