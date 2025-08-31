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

   .. py:method:: get_structure(cell: moyopy.Cell, unique_species_mapping: dict[int, pymatgen.core.Composition] | None = None) -> pymatgen.core.Structure
      :staticmethod:


      Convert a Moyo Cell to a pymatgen Structure.

      If the Cell was created from a disordered Structure, the unique_species_mapping should be
      provided to reconstruct the original species.

      :param cell: The Moyo Cell to convert.
      :type cell: moyopy.Cell
      :param unique_species_mapping: A mapping from integer indices used in the Moyo Cell to the original pymatgen
                                     SpeciesLike objects. If None, assumes the Cell was created from an ordered Structure.
      :type unique_species_mapping: dict[int, Composition] | None

      :returns: **structure** -- The converted pymatgen Structure object.
      :rtype: Structure



   .. py:method:: get_atoms(cell: moyopy.Cell) -> ase.Atoms
      :staticmethod:



   .. py:method:: from_structure(structure: pymatgen.core.Structure) -> moyopy.Cell
      :staticmethod:



   .. py:method:: from_disordered_structure(structure: pymatgen.core.Structure) -> tuple[moyopy.Cell, dict[int, pymatgen.core.Composition]]
      :staticmethod:


      Convert a disordered pymatgen Structure to a Moyo Cell.

      :param structure: A pymatgen Structure object, which may contain disordered sites.
      :type structure: Structure

      :returns: * **cell** (*moyopy.Cell*) -- The converted Moyo Cell object.
                * **unique_species_mapping** (*dict[int, Composition]*) -- A mapping from integer indices used in the Moyo Cell to the original pymatgen
                  SpeciesLike objects.



   .. py:method:: from_atoms(atoms: ase.Atoms) -> moyopy.Cell
      :staticmethod:



   .. py:method:: from_py_obj(struct: pymatgen.core.Structure | ase.Atoms | pymatgen.io.ase.MSONAtoms) -> moyopy.Cell
      :staticmethod:


      Convert a Python atomic structure object to a Moyo Cell.

      :param struct: Currently supports pymatgen Structure, ASE Atoms, and MSONAtoms

      :returns: The converted Moyo cell
      :rtype: moyopy.Cell



