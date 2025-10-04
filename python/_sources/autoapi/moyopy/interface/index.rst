moyopy.interface
================

.. py:module:: moyopy.interface


Classes
-------

.. autoapisummary::

   moyopy.interface.MoyoAdapter
   moyopy.interface.MoyoNonCollinearMagneticAdapter


Functions
---------

.. autoapisummary::

   moyopy.interface._species_from_numbers
   moyopy.interface._get_unique_species_and_mapping


Module Contents
---------------

.. py:class:: MoyoAdapter

   .. py:method:: get_structure(cell: moyopy.Cell, *, unique_species_mapping: dict[int, pymatgen.core.Composition] | None = None) -> pymatgen.core.Structure
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


      Convert a pymatgen Structure to a Moyo Cell.



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



.. py:class:: MoyoNonCollinearMagneticAdapter

   .. py:method:: get_structure(magnetic_cell: moyopy.NonCollinearMagneticCell, *, unique_species_mapping: dict[int, pymatgen.core.Composition] | None = None) -> pymatgen.core.Structure
      :staticmethod:


      Convert a Moyo NonCollinearMagneticCell to a pymatgen Structure.

      If the NonCollinearMagneticCell was created from a disordered Structure, the
      unique_species_mapping should be provided to reconstruct the original species.

      :param magnetic_cell: The Moyo NonCollinearMagneticCell to convert.
      :type magnetic_cell: moyopy.NonCollinearMagneticCell
      :param unique_species_mapping: A mapping from integer indices used in the Moyo NonCollinearMagneticCell to the
                                     original pymatgen SpeciesLike objects. If None, assumes the NonCollinearMagneticCell
                                     was created from an ordered Structure.
      :type unique_species_mapping: dict[int, Composition] | None

      :returns: **structure** -- The converted pymatgen Structure object with magnetic moments in
                site_properties['magmom'].
      :rtype: Structure



   .. py:method:: from_structure(structure: pymatgen.core.Structure) -> moyopy.NonCollinearMagneticCell
      :staticmethod:


      Convert a pymatgen Structure with non-collinear magnetic moments to a Moyo
      NonCollinearMagneticCell.



   .. py:method:: from_disordered_structure(structure: pymatgen.core.Structure) -> tuple[moyopy.NonCollinearMagneticCell, dict[int, pymatgen.core.Composition]]
      :staticmethod:


      Convert a disordered pymatgen Structure with non-collinear magnetic moments to a Moyo
      NonCollinearMagneticCell.

      :param structure: A pymatgen Structure object, which may contain disordered sites. Must have
                        non-collinear magnetic moments in site_properties['magmom'].
      :type structure: Structure

      :returns: * **magnetic_cell** (*moyopy.NonCollinearMagneticCell*) -- The converted Moyo NonCollinearMagneticCell object.
                * **unique_species_mapping** (*dict[int, Composition]*) -- A mapping from integer indices used in the Moyo NonCollinearMagneticCell to the
                  original pymatgen SpeciesLike objects.



.. py:function:: _species_from_numbers(numbers: list[int], unique_species_mapping: dict[int, pymatgen.core.Composition] | None = None) -> list[pymatgen.core.Composition]

.. py:function:: _get_unique_species_and_mapping(structure: pymatgen.core.Structure) -> tuple[dict[pymatgen.core.Composition, int], dict[int, pymatgen.core.Composition]]

