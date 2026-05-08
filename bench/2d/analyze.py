"""Categorize moyopy `MoyoLayerDataset` failures on JARVIS dft_2d.

Companion to ``layer.py``. The benchmark itself swallows exceptions (one
``Abort: jid=... symprec=...`` line per drop) so it can finish the sweep
quickly. This script re-runs the same per-(material, symprec) loop with
full error capture, buckets failures by exception kind + message, and
joins each bucket against the bulk-parent space-group label JARVIS
provides so we can spot systematic groups (e.g. "all bulk-SG-164 failures
have the same root cause").

Outputs:

* A printed report broken down by error-kind, per-symprec totals, and
  per-jid failure pattern (always-fail / partial / always-OK).
* ``stats_moyopy_errors.json`` next to the other bench artifacts, with
  one entry per (jid, symprec) failure carrying the exception message
  and the bulk-parent space group, for follow-up drilling.

Subcommands:

  * ``run``      -- full re-run of the dataset + write JSON.
  * ``report``   -- read an existing JSON and print the summary only.
  * ``probe``    -- drill into one jid: print spglib + moyopy answers at
                    every symprec and the per-symprec error message.

Single material smoke check::

    python analyze.py probe JVASP-76516

Full sweep::

    python analyze.py run                 # writes stats_moyopy_errors.json
    python analyze.py report              # re-print summary from the JSON
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import click
from jarvis.core.atoms import Atoms as JarvisAtoms
from jarvis.db.figshare import data as jarvis_data
from spglib import get_symmetry_layerdataset
from tqdm.auto import tqdm

import moyopy
from moyopy import LayerSetting
from moyopy.interface import MoyoAdapter

SYMPREC_LIST = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
APERIODIC_DIR = 2
DEFAULT_OUT = Path(__file__).resolve().parent / "stats_moyopy_errors.json"


@click.group()
def cli() -> None:
    """Error-classification companion for ``layer.py``."""


@cli.command()
@click.option(
    "--out",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_OUT,
    show_default=True,
    help="Where to write the per-(jid, symprec) failure JSON.",
)
@click.option(
    "--setting",
    type=click.Choice(["spglib", "standard"]),
    default="spglib",
    show_default=True,
    help="LayerSetting to use for MoyoLayerDataset.",
)
def run(out: Path, setting: str) -> None:
    """Re-run the JARVIS dft_2d sweep with full error capture."""
    layer_setting = (
        LayerSetting.spglib() if setting == "spglib" else LayerSetting.standard()
    )
    entries = jarvis_data("dft_2d")

    cell_build_fails: list[dict] = []
    failures: list[dict] = []
    successes = 0
    per_jid_status: dict[str, dict[float, str]] = {}

    for entry in tqdm(entries):
        jid = entry["jid"]
        spg_num = entry.get("spg_number", "?")
        formula = entry.get("formula", "?")
        atoms = JarvisAtoms.from_dict(entry["atoms"])
        structure = atoms.pymatgen_converter()
        try:
            cell = MoyoAdapter.from_structure(structure)
        except Exception as e:
            cell_build_fails.append(
                {
                    "jid": jid,
                    "exception": type(e).__name__,
                    "msg": str(e).strip() or "<empty>",
                    "bulk_spg": spg_num,
                    "formula": formula,
                }
            )
            continue

        per_jid_status[jid] = {}
        for symprec in SYMPREC_LIST:
            try:
                moyopy.MoyoLayerDataset(cell, symprec=symprec, setting=layer_setting)
                successes += 1
                per_jid_status[jid][symprec] = "OK"
            except Exception as e:
                kind = type(e).__name__
                msg = str(e).strip() or "<empty>"
                key = f"{kind}: {msg}"
                failures.append(
                    {
                        "jid": jid,
                        "symprec": symprec,
                        "exception": kind,
                        "msg": msg,
                        "bulk_spg": spg_num,
                        "formula": formula,
                    }
                )
                per_jid_status[jid][symprec] = key

    payload = {
        "n_entries": len(entries),
        "n_symprec": len(SYMPREC_LIST),
        "symprec_list": SYMPREC_LIST,
        "setting": setting,
        "successes": successes,
        "cell_build_fails": cell_build_fails,
        "failures": failures,
    }
    out.write_text(json.dumps(payload, indent=2))
    click.echo(f"Wrote {out}")
    _print_report(payload)


@cli.command()
@click.option(
    "--src",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_OUT,
    show_default=True,
)
def report(src: Path) -> None:
    """Re-print the summary from an existing JSON without re-running."""
    payload = json.loads(src.read_text())
    _print_report(payload)


@cli.command()
@click.argument("jid")
@click.option(
    "--setting",
    type=click.Choice(["spglib", "standard"]),
    default="spglib",
    show_default=True,
)
def probe(jid: str, setting: str) -> None:
    """Drill into a single material: spglib vs moyopy at every symprec."""
    layer_setting = (
        LayerSetting.spglib() if setting == "spglib" else LayerSetting.standard()
    )
    entries = jarvis_data("dft_2d")
    entry = next((e for e in entries if e["jid"] == jid), None)
    if entry is None:
        raise click.ClickException(f"jid {jid!r} not found in JARVIS dft_2d")
    atoms = JarvisAtoms.from_dict(entry["atoms"])
    click.echo(f"jid      : {jid}")
    click.echo(f"formula  : {entry['formula']}")
    click.echo(f"nat      : {entry['nat']}")
    click.echo(f"bulk SG  : {entry['spg_number']} {entry['spg_symbol']}")
    click.echo(f"crys     : {entry.get('crys', '?')}")
    click.echo("lattice (rows = basis vectors):")
    for row in atoms.lattice_mat:
        click.echo(f"  {[float(x) for x in row]}")
    click.echo("frac_coords:")
    for sym, p in zip(atoms.elements, atoms.frac_coords):
        click.echo(f"  {sym}  {[float(x) for x in p]}")

    spg_cell = (atoms.lattice_mat, atoms.frac_coords, atoms.atomic_numbers)
    click.echo("\nspglib:")
    for sp in SYMPREC_LIST:
        try:
            ds = get_symmetry_layerdataset(spg_cell, aperiodic_dir=APERIODIC_DIR, symprec=sp)
            click.echo(
                f"  symprec={sp}: LG={ds.number} hall={ds.hall_number} "
                f"n_ops={len(ds.rotations)}"
            )
        except Exception as e:
            click.echo(f"  symprec={sp}: ERR {type(e).__name__}: {e}")

    structure = atoms.pymatgen_converter()
    try:
        cell = MoyoAdapter.from_structure(structure)
    except Exception as e:
        click.echo(f"\nmoyopy cell build failed: {type(e).__name__}: {e}")
        return
    click.echo(f"\nmoyopy (setting={setting}):")
    for sp in SYMPREC_LIST:
        try:
            ds = moyopy.MoyoLayerDataset(cell, symprec=sp, setting=layer_setting)
            click.echo(
                f"  symprec={sp}: LG={ds.number} hall={ds.hall_number} "
                f"n_ops={len(ds.operations.rotations)}"
            )
        except Exception as e:
            click.echo(f"  symprec={sp}: ERR {type(e).__name__}: {e}")


def _print_report(payload: dict) -> None:
    failures = payload["failures"]
    cell_build_fails = payload["cell_build_fails"]
    n_entries = payload["n_entries"]
    n_symprec = payload["n_symprec"]
    successes = payload["successes"]
    total_attempts = n_entries * n_symprec - len(cell_build_fails) * n_symprec

    click.secho(
        f"\n=== Cell-build failures (MoyoAdapter.from_structure): {len(cell_build_fails)} ===",
        bold=True,
    )
    if cell_build_fails:
        c = Counter(f["exception"] for f in cell_build_fails)
        for kind, n in c.most_common():
            click.echo(f"  {kind}: {n}")

    click.secho(
        f"\n=== Compute failures: {len(failures)} / {total_attempts} "
        f"({len(failures) / total_attempts * 100:.1f}%); successes: {successes} ===",
        bold=True,
    )

    by_kind: dict[str, list[dict]] = defaultdict(list)
    for f in failures:
        by_kind[f"{f['exception']}: {f['msg']}"].append(f)
    ranked = sorted(by_kind.items(), key=lambda kv: -len(kv[1]))

    click.echo("\nTop 15 error kinds (by total count):")
    for key, occs in ranked[:15]:
        n = len(occs)
        n_jids = len({f["jid"] for f in occs})
        bulk_spg_top = Counter(f["bulk_spg"] for f in occs).most_common(3)
        bulk_str = ", ".join(f"SG{sp}={k}" for sp, k in bulk_spg_top)
        click.echo(
            f"  [{n:>4}] in {n_jids:>3} unique jids  bulk-SG: {bulk_str}"
        )
        click.echo(f"          msg: {key}")
        for sample in occs[:2]:
            click.echo(
                f"          ex: {sample['jid']} symprec={sample['symprec']} "
                f"bulk_spg={sample['bulk_spg']} {sample['formula']}"
            )

    click.echo("\nPer-symprec failure counts:")
    by_sp: dict[float, int] = defaultdict(int)
    for f in failures:
        by_sp[f["symprec"]] += 1
    for sp in payload["symprec_list"]:
        click.echo(f"  symprec={sp}: {by_sp.get(sp, 0)} failures")

    # Per-jid pattern: count materials by always-fail / partial / always-OK.
    by_jid: dict[str, set[float]] = defaultdict(set)
    all_jids: set[str] = set()
    for f in failures:
        by_jid[f["jid"]].add(f["symprec"])
        all_jids.add(f["jid"])
    n_symprec_set = set(payload["symprec_list"])
    full_fail = sum(1 for j, sps in by_jid.items() if sps == n_symprec_set)
    partial_fail = sum(1 for j, sps in by_jid.items() if sps != n_symprec_set)

    n_seen = n_entries - len(cell_build_fails)
    click.echo("\nPer-jid failure pattern:")
    click.echo(f"  {n_seen - len(by_jid)} jids: OK at every symprec")
    click.echo(f"  {partial_fail} jids: OK at some, fail at others")
    click.echo(f"  {full_fail} jids: fail at every symprec")

    # Per-error-kind, per-bulk-SG cross-tab (top kind only).
    if ranked:
        top_key, top_occs = ranked[0]
        sg_counts = Counter(f["bulk_spg"] for f in top_occs)
        click.secho(f"\nBulk-SG breakdown of dominant error ({top_key}):", bold=True)
        for sg, n in sg_counts.most_common(10):
            n_jids = len({f["jid"] for f in top_occs if f["bulk_spg"] == sg})
            click.echo(f"  bulk SG {sg}: {n} (in {n_jids} jids)")


def _iter_failure_jids(failures: Iterable[dict]) -> Iterable[str]:
    seen: set[str] = set()
    for f in failures:
        if f["jid"] not in seen:
            seen.add(f["jid"])
            yield f["jid"]


if __name__ == "__main__":
    cli()
