from __future__ import annotations

import json
from time import perf_counter

import matbench_discovery.data
from pymatgen.entries.computed_entries import ComputedStructureEntry
from spglib import get_symmetry_dataset
from tqdm.auto import tqdm

import moyopy


def main():
    df = matbench_discovery.data.load("mp_computed_structure_entries", version="1.0.0")

    symprec = 1e-4

    time_moyopy = 0
    time_spglib = 0
    with tqdm(df.iterrows(), total=len(df)) as pbar:
        for material_id, row in pbar:
            structure = ComputedStructureEntry.from_dict(row["entry"]).structure
            basis = structure.lattice.matrix
            positions = structure.frac_coords
            numbers = [site.specie.Z for site in structure]

            try:
                moyopy_cell = moyopy.Cell(basis.tolist(), positions.tolist(), numbers)
                start = perf_counter()
                moyopy_dataset = moyopy.MoyoDataset(moyopy_cell, symprec=symprec)
                time_moyopy += perf_counter() - start

                spglib_cell = (basis, positions, numbers)
                start = perf_counter()
                spglib_dataset = get_symmetry_dataset(spglib_cell, symprec=symprec)
                time_spglib += perf_counter() - start
                assert moyopy_dataset.number == spglib_dataset["number"]
                pbar.set_postfix_str(f"{material_id=} number={moyopy_dataset.number}")
            except:  # noqa: E722
                print(f"Failed: {material_id=}")
                print(f"Elapsed time: {time_moyopy=}, {time_spglib=}")

                with open(f"{material_id}.json", "w") as f:
                    json.dump(structure.as_dict(), f, indent=4)
                break


if __name__ == "__main__":
    main()
