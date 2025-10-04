from __future__ import annotations

from warnings import warn

try:
    import ase
    import numpy as np
    from pymatgen.core import Composition, Element, Structure
    from pymatgen.io.ase import MSONAtoms
except ImportError:
    raise ImportError("Try installing dependencies with `pip install moyopy[interface]`")

import moyopy


class MoyoAdapter:
    @staticmethod
    def get_structure(
        cell: moyopy.Cell, *, unique_species_mapping: dict[int, Composition] | None = None
    ) -> Structure:
        """Convert a Moyo Cell to a pymatgen Structure.

        If the Cell was created from a disordered Structure, the unique_species_mapping should be
        provided to reconstruct the original species.

        Parameters
        ----------
        cell : moyopy.Cell
            The Moyo Cell to convert.
        unique_species_mapping : dict[int, Composition] | None
            A mapping from integer indices used in the Moyo Cell to the original pymatgen
            SpeciesLike objects. If None, assumes the Cell was created from an ordered Structure.

        Returns
        -------
        structure : Structure
            The converted pymatgen Structure object.
        """
        species = _species_from_numbers(cell.numbers, unique_species_mapping)
        return Structure(lattice=cell.basis, species=species, coords=cell.positions)

    @staticmethod
    def get_atoms(cell: moyopy.Cell) -> ase.Atoms:
        atoms = ase.Atoms(
            cell=cell.basis,
            scaled_positions=cell.positions,
            numbers=cell.numbers,
            pbc=True,
        )
        return atoms

    @staticmethod
    def from_structure(structure: Structure) -> moyopy.Cell:
        """Convert a pymatgen Structure to a Moyo Cell."""
        if not structure.is_ordered:
            raise ValueError("Structure must be ordered. Use from_disordered_structure instead.")

        basis = structure.lattice.matrix.tolist()
        positions = structure.frac_coords.tolist()
        numbers = [site.specie.Z for site in structure]

        return moyopy.Cell(
            basis=basis,
            positions=positions,
            numbers=numbers,
        )

    @staticmethod
    def from_disordered_structure(
        structure: Structure,
    ) -> tuple[moyopy.Cell, dict[int, Composition]]:
        """Convert a disordered pymatgen Structure to a Moyo Cell.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object, which may contain disordered sites.

        Returns
        -------
        cell : moyopy.Cell
            The converted Moyo Cell object.
        unique_species_mapping : dict[int, Composition]
            A mapping from integer indices used in the Moyo Cell to the original pymatgen
            SpeciesLike objects.
        """
        if structure.is_ordered:
            warn("Structure is ordered. Consider using from_structure.")

        unique_species, unique_species_mapping = _get_unique_species_and_mapping(structure)

        basis = structure.lattice.matrix.tolist()
        positions = structure.frac_coords.tolist()
        numbers = [unique_species[site.species] for site in structure]
        cell = moyopy.Cell(
            basis=basis,
            positions=positions,
            numbers=numbers,
        )

        return cell, unique_species_mapping

    @staticmethod
    def from_atoms(atoms: ase.Atoms) -> moyopy.Cell:
        basis = list(atoms.cell)
        positions = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()

        return moyopy.Cell(
            basis=basis,
            positions=positions,
            numbers=numbers,
        )

    @staticmethod
    def from_py_obj(struct: Structure | ase.Atoms | MSONAtoms) -> moyopy.Cell:
        """Convert a Python atomic structure object to a Moyo Cell.

        Args:
            struct: Currently supports pymatgen Structure, ASE Atoms, and MSONAtoms

        Returns:
            moyopy.Cell: The converted Moyo cell
        """
        if isinstance(struct, (ase.Atoms, MSONAtoms)):
            return MoyoAdapter.from_atoms(struct)
        elif isinstance(struct, Structure):
            return MoyoAdapter.from_structure(struct)
        else:
            cls_name = type(struct).__name__
            raise TypeError(f"Expected Structure, Atoms, or MSONAtoms, got {cls_name}")


class MoyoNonCollinearMagneticAdapter:
    @staticmethod
    def get_structure(
        magnetic_cell: moyopy.NonCollinearMagneticCell,
        *,
        unique_species_mapping: dict[int, Composition] | None = None,
    ) -> Structure:
        """Convert a Moyo NonCollinearMagneticCell to a pymatgen Structure.

        If the NonCollinearMagneticCell was created from a disordered Structure, the
        unique_species_mapping should be provided to reconstruct the original species.

        Parameters
        ----------
        magnetic_cell : moyopy.NonCollinearMagneticCell
            The Moyo NonCollinearMagneticCell to convert.
        unique_species_mapping : dict[int, Composition] | None
            A mapping from integer indices used in the Moyo NonCollinearMagneticCell to the
            original pymatgen SpeciesLike objects. If None, assumes the NonCollinearMagneticCell
            was created from an ordered Structure.

        Returns
        -------
        structure : Structure
            The converted pymatgen Structure object with magnetic moments in
            site_properties['magmom'].
        """
        species = _species_from_numbers(magnetic_cell.numbers, unique_species_mapping)
        return Structure(
            lattice=magnetic_cell.basis,
            species=species,
            coords=magnetic_cell.positions,
            site_properties={"magmom": magnetic_cell.magnetic_moments},
        )

    @staticmethod
    def from_structure(structure: Structure) -> moyopy.NonCollinearMagneticCell:
        """Convert a pymatgen Structure with non-collinear magnetic moments to a Moyo
        NonCollinearMagneticCell."""
        if not structure.is_ordered:
            raise ValueError("Structure must be ordered. Use from_disordered_structure instead.")

        basis = structure.lattice.matrix.tolist()
        positions = structure.frac_coords.tolist()
        numbers = [site.specie.Z for site in structure]
        if "magmom" not in structure.site_properties:
            raise ValueError(
                "Structure must have non-collinear magnetic moments in site_properties['magmom']"
            )
        magnetic_moments_arr = np.array(structure.site_properties["magmom"])

        if magnetic_moments_arr.shape != (len(structure), 3):
            raise ValueError(
                f"Structure must have non-collinear magnetic moments. "
                f"Expected: ({len(structure)}, 3), actual: {magnetic_moments_arr.shape}."
            )

        return moyopy.NonCollinearMagneticCell(
            basis=basis,
            positions=positions,
            numbers=numbers,
            magnetic_moments=magnetic_moments_arr.tolist(),
        )

    @staticmethod
    def from_disordered_structure(
        structure: Structure,
    ) -> tuple[moyopy.NonCollinearMagneticCell, dict[int, Composition]]:
        """Convert a disordered pymatgen Structure with non-collinear magnetic moments to a Moyo
        NonCollinearMagneticCell.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object, which may contain disordered sites. Must have
            non-collinear magnetic moments in site_properties['magmom'].

        Returns
        -------
        magnetic_cell : moyopy.NonCollinearMagneticCell
            The converted Moyo NonCollinearMagneticCell object.
        unique_species_mapping : dict[int, Composition]
            A mapping from integer indices used in the Moyo NonCollinearMagneticCell to the
            original pymatgen SpeciesLike objects.
        """
        if structure.is_ordered:
            warn("Structure is ordered. Consider using from_structure.")

        unique_species, unique_species_mapping = _get_unique_species_and_mapping(structure)

        basis = structure.lattice.matrix.tolist()
        positions = structure.frac_coords.tolist()
        numbers = [unique_species[site.species] for site in structure]
        if "magmom" not in structure.site_properties:
            raise ValueError(
                "Structure must have non-collinear magnetic moments in site_properties['magmom']"
            )
        magnetic_moments_arr = np.array(structure.site_properties["magmom"])

        if magnetic_moments_arr.shape != (len(structure), 3):
            raise ValueError(
                "Structure must have non-collinear magnetic moments in site_properties['magmom']. "
                f"Expected shape: ({len(structure)}, 3), got: {magnetic_moments_arr.shape}."
            )

        magnetic_cell = moyopy.NonCollinearMagneticCell(
            basis=basis,
            positions=positions,
            numbers=numbers,
            magnetic_moments=magnetic_moments_arr.tolist(),
        )
        return magnetic_cell, unique_species_mapping


def _species_from_numbers(
    numbers: list[int], unique_species_mapping: dict[int, Composition] | None = None
) -> list[Composition]:
    if unique_species_mapping:
        species = [unique_species_mapping[number] for number in numbers]
    else:
        species = [Element.from_Z(number) for number in numbers]
    return species


def _get_unique_species_and_mapping(
    structure: Structure,
) -> tuple[dict[Composition, int], dict[int, Composition]]:
    unique_species = {}
    for site in structure:
        key = site.species
        if key not in unique_species:
            unique_species[key] = len(unique_species)
    unique_species_mapping = {i: species for i, species in enumerate(unique_species)}
    return unique_species, unique_species_mapping
