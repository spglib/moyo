from __future__ import annotations

from time import perf_counter

import click
import pandas as pd
from jarvis.core.atoms import Atoms as JarvisAtoms
from jarvis.db.figshare import data as jarvis_data
from spglib import get_symmetry_layerdataset
from tqdm.auto import tqdm

import moyopy
from moyopy import LayerSetting
from moyopy.interface import MoyoAdapter

SYMPREC_LIST = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
APERIODIC_DIR = 2


@click.group()
def cli(): ...


@cli.command()
def perf_spglib():
    entries = jarvis_data("dft_2d")
    all_stats = []

    with tqdm(entries) as pbar:
        for entry in pbar:
            jid = entry["jid"]
            atoms = JarvisAtoms.from_dict(entry["atoms"])
            basis = atoms.lattice_mat
            positions = atoms.frac_coords
            numbers = atoms.atomic_numbers
            spglib_cell = (basis, positions, numbers)

            for i, symprec in enumerate(SYMPREC_LIST):
                try:
                    start = perf_counter()
                    spglib_dataset = get_symmetry_layerdataset(
                        spglib_cell,
                        aperiodic_dir=APERIODIC_DIR,
                        symprec=symprec,
                    )
                    time_spglib = perf_counter() - start

                    all_stats.append(
                        {
                            "id": f"{jid}_{i}",
                            "jid": jid,
                            "symprec": symprec,
                            "time_spglib": time_spglib,
                            "num_atoms": len(numbers),
                            "number_spglib": spglib_dataset["number"],
                        }
                    )
                except:  # noqa: E722
                    print(f"Abort: {jid=} {symprec=}")

    df_stats = pd.DataFrame(all_stats)
    df_stats.to_json("stats_spglib.json")


@cli.command()
def perf_moyopy():
    entries = jarvis_data("dft_2d")
    all_stats = []

    with tqdm(entries) as pbar:
        for entry in pbar:
            jid = entry["jid"]
            atoms = JarvisAtoms.from_dict(entry["atoms"])
            structure = atoms.pymatgen_converter()
            try:
                moyopy_cell = MoyoAdapter.from_structure(structure)
            except:  # noqa: E722
                print(f"Skip (cell build failed): {jid=}")
                continue

            for i, symprec in enumerate(SYMPREC_LIST):
                try:
                    start = perf_counter()
                    moyopy_dataset = moyopy.MoyoLayerDataset(
                        moyopy_cell,
                        symprec=symprec,
                        setting=LayerSetting.spglib(),
                    )
                    time_moyopy = perf_counter() - start

                    all_stats.append(
                        {
                            "id": f"{jid}_{i}",
                            "jid": jid,
                            "symprec": symprec,
                            "time_moyopy": time_moyopy,
                            "num_atoms": len(structure),
                            "number_moyopy": moyopy_dataset.number,
                        }
                    )
                except:  # noqa: E722
                    print(f"Abort: {jid=} {symprec=}")

    df_stats = pd.DataFrame(all_stats)
    df_stats.to_json("stats_moyopy.json")


if __name__ == "__main__":
    cli()
