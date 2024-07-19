from __future__ import annotations

try:
    import ase
    from pymatgen.core import Element, Structure
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
