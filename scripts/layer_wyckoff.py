"""Generate moyo's layer-group Wyckoff position database.

Fetches `database/layer_Wyckoff.csv` from spglib's GitHub (default: tip of
`develop`; pin with `--ref <SHA>` for reproducible regeneration) and emits
Rust source for the `LayerWyckoffPosition` array.

Usage:
    python scripts/layer_wyckoff.py [--ref <git ref or sha>]

The CSV is keyed by spglib's layer Hall number (1-116), the same numbering
that `database/layer_spg.csv` uses. moyo's `LayerHallNumber` is row-aligned
with that file, so the Hall numbers carry over directly.

CSV format (one Wyckoff position per logical row):
    `<hall>:<HM>:::::::`                                           -- LG header
    `::<mult>:<letter>:<sym>:<coord1>:<coord2>:<coord3>:<coord4>`  -- position
    `:::::<coord5>:<coord6>:<coord7>:<coord8>`                     -- continuation

We extract only the representative coordinate (`coord1`) per position;
`assign_wyckoff_position` matches Wyckoff orbits via integer offset search,
so we don't need the rest.
"""

import io
import sys
import urllib.request

import click

REPO = "spglib/spglib"
CSV_PATH = "database/layer_Wyckoff.csv"
DEFAULT_REF = "develop"


def fetch_csv(ref: str) -> str:
    url = f"https://raw.githubusercontent.com/{REPO}/{ref}/{CSV_PATH}"
    with urllib.request.urlopen(url) as resp:
        return resp.read().decode("utf-8")


@click.command()
@click.option(
    "--ref",
    default=DEFAULT_REF,
    show_default=True,
    help="git ref or full SHA to fetch from spglib.",
)
def main(ref: str) -> None:
    raw = fetch_csv(ref)

    rows = []
    hall_number = None
    for line in io.StringIO(raw):
        line = line.strip("\n")
        if not line:
            continue
        if line == "end of data":
            break
        cols = line.split(":")
        # Header row: `<hall>:<HM>:::::::` -- column 0 nonempty (hall number).
        if cols[0] != "":
            hall_number = int(cols[0])
            continue
        # Continuation row: `:::::<coord5>:<coord6>:...` -- column 2 empty.
        if cols[2] == "":
            continue
        # Position row: `::<mult>:<letter>:<sym>:<coord1>:<coord2>:<coord3>:<coord4>`.
        multiplicity = int(cols[2])
        letter = cols[3]
        site_symmetry = cols[4]
        coord1 = cols[5].strip("(").strip(")")
        rows.append((hall_number, multiplicity, letter, site_symmetry, coord1))

    print(
        f"// Source: https://github.com/{REPO}/blob/{ref}/{CSV_PATH}",
        file=sys.stderr,
    )
    print(f"// {len(rows)} entries", file=sys.stderr)

    print(
        f"const LAYER_WYCKOFF_DATABASE: [LayerWyckoffPosition; {len(rows)}] = ["
    )
    for hall_number, multiplicity, letter, site_symmetry, coord in rows:
        print(
            f"    LayerWyckoffPosition::new({hall_number}, {multiplicity}, "
            f"'{letter}', \"{site_symmetry}\", \"{coord}\"),"
        )
    print("];")


if __name__ == "__main__":
    main()
