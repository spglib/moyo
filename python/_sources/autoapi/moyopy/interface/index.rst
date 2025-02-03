moyopy.interface
================

.. py:module:: moyopy.interface


Classes
-------

.. autoapisummary::

   moyopy.interface.MoyoAdapter


Module Contents
---------------

.. py:class:: MoyoAdapter

   .. py:method:: get_structure(cell: moyopy.Cell) -> pymatgen.core.Structure
      :staticmethod:



   .. py:method:: get_atoms(cell: moyopy.Cell) -> ase.Atoms
      :staticmethod:



   .. py:method:: from_structure(structure: pymatgen.core.Structure) -> moyopy.Cell
      :staticmethod:



   .. py:method:: from_atoms(atoms: ase.Atoms) -> moyopy.Cell
      :staticmethod:



   .. py:method:: from_py_obj(struct: pymatgen.core.Structure | ase.Atoms | pymatgen.io.ase.MSONAtoms) -> moyopy.Cell
      :staticmethod:


      Convert a Python atomic structure object to a Moyo Cell.

      :param struct: Currently supports pymatgen Structure, ASE Atoms, and MSONAtoms

      :returns: The converted Moyo cell
      :rtype: moyopy.Cell



