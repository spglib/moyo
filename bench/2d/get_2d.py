from __future__ import annotations

from jarvis.core.atoms import Atoms as JarvisAtoms
from jarvis.db.figshare import data as jarvis_data
from spglib import get_symmetry_layerdataset

import moyopy
from moyopy.interface import MoyoAdapter


def main():
    entries = jarvis_data("dft_2d")
    sample = entries[0]
    jid = sample["jid"]
    atoms = JarvisAtoms.from_dict(sample["atoms"])

    structure = atoms.pymatgen_converter()
    cell = MoyoAdapter.from_structure(structure)
    spg_cell = (atoms.lattice_mat, atoms.frac_coords, atoms.atomic_numbers)

    print(f"{jid=} num_atoms={len(structure)}")
    for symprec in [1e-4, 1e-3, 1e-2, 1e-1]:
        ds_s = get_symmetry_layerdataset(spg_cell, aperiodic_dir=2, symprec=symprec)
        ds_m = moyopy.MoyoLayerDataset(cell, symprec=symprec)
        print(symprec, ds_s["number"], ds_m.number)


if __name__ == "__main__":
    main()
