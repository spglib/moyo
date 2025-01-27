from __future__ import annotations

from time import perf_counter

import click
import pandas as pd
from matbench_discovery.data import DataFiles
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key
from spglib import get_symmetry_dataset
from tqdm.auto import tqdm

import moyopy
from moyopy.interface import MoyoAdapter

SYMPREC_LIST = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]


@click.group()
def cli(): ...


@cli.command()
def perf_spglib():
    data_path = DataFiles.mp_computed_structure_entries.path
    df = pd.read_json(data_path).set_index(Key.mat_id)

    all_stats = []

    with tqdm(df.iterrows(), total=len(df)) as pbar:
        for material_id, row in pbar:
            structure = ComputedStructureEntry.from_dict(row["entry"]).structure
            basis = structure.lattice.matrix
            positions = structure.frac_coords
            numbers = [site.specie.Z for site in structure]
            spglib_cell = (basis, positions, numbers)

            for i, symprec in enumerate(SYMPREC_LIST):
                try:
                    start = perf_counter()
                    spglib_dataset = get_symmetry_dataset(spglib_cell, symprec=symprec)
                    time_spglib = perf_counter() - start

                    all_stats.append(
                        {
                            "id": f"{material_id}_{i}",
                            "material_id": material_id,
                            "symprec": symprec,
                            "time_spglib": time_spglib,
                            "num_atoms": len(numbers),
                            "number_spglib": spglib_dataset["number"],
                        }
                    )
                except:  # noqa: E722
                    print(f"Abort: {material_id=} {symprec=}")

    df_stats = pd.DataFrame(all_stats)
    df_stats.to_json("stats_spglib.json")


@cli.command()
def perf_moyopy():
    data_path = DataFiles.mp_computed_structure_entries.path
    df = pd.read_json(data_path).set_index(Key.mat_id)

    all_stats = []

    with tqdm(df.iterrows(), total=len(df)) as pbar:
        for material_id, row in pbar:
            structure = ComputedStructureEntry.from_dict(row["entry"]).structure
            numbers = [site.specie.Z for site in structure]
            moyopy_cell = MoyoAdapter.from_structure(structure)

            for i, symprec in enumerate(SYMPREC_LIST):
                try:
                    start = perf_counter()
                    moyopy_dataset = moyopy.MoyoDataset(moyopy_cell, symprec=symprec)
                    time_moyopy = perf_counter() - start

                    all_stats.append(
                        {
                            "id": f"{material_id}_{i}",
                            "material_id": material_id,
                            "symprec": symprec,
                            "time_moyopy": time_moyopy,
                            "num_atoms": len(numbers),
                            "number_moyopy": moyopy_dataset.number,
                        }
                    )
                except:  # noqa: E722
                    print(f"Abort: {material_id=} {symprec=}")

                    with open(f"{material_id}.json", "w") as f:
                        f.write(moyopy_cell.serialize_json())
                    return

    df_stats = pd.DataFrame(all_stats)
    df_stats.to_json("stats_moyopy.json")


if __name__ == "__main__":
    cli()
