from __future__ import annotations

try:
    import ase
    from pymatgen.core import Element, Structure
    from pymatgen.io.ase import MSONAtoms
except ImportError:
    raise ImportError("Try installing dependencies with `pip install moyopy[interface]`")

import moyopy


class MoyoAdapter:
    @staticmethod
    def get_structure(cell: moyopy.Cell) -> Structure:
        species = [Element.from_Z(number) for number in cell.numbers]
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
        basis = structure.lattice.matrix
        positions = structure.frac_coords
        numbers = [site.specie.Z for site in structure]

        return moyopy.Cell(
            basis=basis.tolist(),
            positions=positions.tolist(),
            numbers=numbers,
        )

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
