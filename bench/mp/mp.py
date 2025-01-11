from __future__ import annotations

from time import perf_counter

import pandas as pd
import spglib
from matbench_discovery.data import DataFiles
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key
from spglib import get_symmetry_dataset
from tqdm.auto import tqdm

import moyopy


def main():
    data_path = DataFiles.mp_computed_structure_entries.path
    df = pd.read_json(data_path).set_index(Key.mat_id)

    all_stats = []

    with tqdm(df.iterrows(), total=len(df)) as pbar:
        for material_id, row in pbar:
            structure = ComputedStructureEntry.from_dict(row["entry"]).structure
            basis = structure.lattice.matrix
            positions = structure.frac_coords
            numbers = [site.specie.Z for site in structure]

            moyopy_cell = moyopy.Cell(basis.tolist(), positions.tolist(), numbers)
            spglib_cell = (basis, positions, numbers)

            for symprec in [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]:
                try:
                    start = perf_counter()
                    moyopy_dataset = moyopy.MoyoDataset(moyopy_cell, symprec=symprec)
                    time_moyopy = perf_counter() - start

                    start = perf_counter()
                    spglib_dataset = get_symmetry_dataset(spglib_cell, symprec=symprec)
                    time_spglib = perf_counter() - start

                    all_stats.append(
                        {
                            "time_moyopy": time_moyopy,
                            "time_spglib": time_spglib,
                            "symprec": symprec,
                            "num_atoms": len(numbers),
                            "material_id": material_id,
                            "number_moyopy": moyopy_dataset.number,
                            "number_spglib": spglib_dataset["number"],
                            "moyopy": moyopy.__version__,
                            "spglib": spglib.__version__,
                        }
                    )
                except:  # noqa: E722
                    print(f"Abort: {material_id=} {symprec=}")

                    with open(f"{material_id}.json", "w") as f:
                        f.write(moyopy_cell.serialize_json())
                    return

                # if moyopy_dataset.number != spglib_dataset["number"]:
                #     print(
                #         f"Inconsistent: {material_id=} {moyopy_dataset.number=} {spglib_dataset['number']=}"  # noqa: E501
                #     )

    df_stats = pd.DataFrame(all_stats)
    df_stats.to_json("stats.json")


if __name__ == "__main__":
    main()
