"""Generate moyo's layer-group Hall symbol database.

Fetches `database/layer_spg.csv` from spglib's GitHub (default: tip of
`develop`; pin with `--ref <SHA>` for reproducible regeneration) and emits
Rust source for the `LayerHallSymbolEntry` array.

Usage:
    python scripts/layer_hall_symbol.py [--ref <git ref or sha>]

Mirrors `scripts/hall_symbol.py` (the 3D space-group analogue) but downloads
its input over HTTP rather than taking a local CSV path. spglib's
`layer_spg.csv` carries one row per Hall setting (~116 rows covering the 80
layer groups). The mapping from layer-group number to layer arithmetic
crystal class is hard-coded here from the boundaries documented in
`moyo/src/data/layer_arithmetic_crystal_class.rs` (paper Fu et al. 2024
Table 4) -- spglib does not publish that mapping in machine-readable form.
"""

import csv
import io
import sys
import urllib.request

import click

REPO = "spglib/spglib"
CSV_PATH = "database/layer_spg.csv"
DEFAULT_REF = "develop"

# Layer arithmetic crystal classes (paper Fu et al. 2024 Table 4): each class
# is a contiguous run of layer-group numbers. Boundaries below come from the
# `// LG <n>` smallest-LG comments in
# `moyo/src/data/layer_arithmetic_crystal_class.rs`; cross-checked against
# the (Schoenflies, lattice symbol) of every row in `layer_spg.csv`.
_ARITH_CLASS_RANGES = [
    (1, 1),  # 1   p1       (C1,  mp)
    (2, 2),  # 2   p-1      (Ci,  mp)
    (3, 3),  # 3   p112     (C2,  mp)
    (4, 5),  # 4   p11m     (C1h, mp)
    (6, 7),  # 5   p112/m   (C2h, mp)
    (8, 9),  # 6   p211     (C2,  op)
    (10, 10),  # 7   c211     (C2,  oc)
    (11, 12),  # 8   pm11     (C1h, op)
    (13, 13),  # 9   cm11     (C1h, oc)
    (14, 17),  # 10  p2/m11   (C2h, op)
    (18, 18),  # 11  c2/m11   (C2h, oc)
    (19, 21),  # 12  p222     (D2,  op)
    (22, 22),  # 13  c222     (D2,  oc)
    (23, 25),  # 14  pmm2     (C2v, op)
    (26, 26),  # 15  cmm2     (C2v, oc)
    (27, 34),  # 16  pm2m     (C2v, op)
    (35, 36),  # 17  cm2m     (C2v, oc)
    (37, 46),  # 18  pmmm     (D2h, op)
    (47, 48),  # 19  cmmm     (D2h, oc)
    (49, 49),  # 20  p4       (C4,  tp)
    (50, 50),  # 21  p-4      (S4,  tp)
    (51, 52),  # 22  p4/m     (C4h, tp)
    (53, 54),  # 23  p422     (D4,  tp)
    (55, 56),  # 24  p4mm     (C4v, tp)
    (57, 58),  # 25  p-42m    (D2d, tp)
    (59, 60),  # 26  p-4m2    (D2d, tp)
    (61, 64),  # 27  p4/mmm   (D4h, tp)
    (65, 65),  # 28  p3       (C3,  hp)
    (66, 66),  # 29  p-3      (C3i, hp)
    (67, 67),  # 30  p312     (D3,  hp)
    (68, 68),  # 31  p321     (D3,  hp)
    (69, 69),  # 32  p3m1     (C3v, hp)
    (70, 70),  # 33  p31m     (C3v, hp)
    (71, 71),  # 34  p-31m    (D3d, hp)
    (72, 72),  # 35  p-3m1    (D3d, hp)
    (73, 73),  # 36  p6       (C6,  hp)
    (74, 74),  # 37  p-6      (C3h, hp)
    (75, 75),  # 38  p6/m     (C6h, hp)
    (76, 76),  # 39  p622     (D6,  hp)
    (77, 77),  # 40  p6mm     (C6v, hp)
    (78, 78),  # 41  p-6m2    (D3h, hp)
    (79, 79),  # 42  p-62m    (D3h, hp)
    (80, 80),  # 43  p6/mmm   (D6h, hp)
]

assert len(_ARITH_CLASS_RANGES) == 43

LG_TO_ARITH: dict[int, int] = {}
for _arith, (_lo, _hi) in enumerate(_ARITH_CLASS_RANGES, start=1):
    for _lg in range(_lo, _hi + 1):
        LG_TO_ARITH[_lg] = _arith
assert sorted(LG_TO_ARITH.keys()) == list(range(1, 81))


def fetch_csv(ref: str) -> str:
    url = f"https://raw.githubusercontent.com/{REPO}/{ref}/{CSV_PATH}"
    with urllib.request.urlopen(url) as resp:
        return resp.read().decode("utf-8")


def centering_for(hall_symbol: str) -> str:
    head = hall_symbol.lstrip("-").lstrip().split(" ", 1)[0]
    if head == "p":
        return "LayerCentering::P"
    if head == "c":
        return "LayerCentering::C"
    raise ValueError(f"unexpected layer Hall lattice symbol: {hall_symbol!r}")


@click.command()
@click.option(
    "--ref",
    default=DEFAULT_REF,
    show_default=True,
    help="git ref or full SHA to fetch from spglib.",
)
def main(ref: str) -> None:
    raw = fetch_csv(ref)
    rows = [row for row in csv.reader(io.StringIO(raw)) if row]

    print(
        f"// Source: https://github.com/{REPO}/blob/{ref}/{CSV_PATH}",
        file=sys.stderr,
    )
    print(f"// {len(rows)} entries", file=sys.stderr)

    print(
        f"const LAYER_HALL_SYMBOL_DATABASE: [LayerHallSymbolEntry; {len(rows)}] = ["
    )
    for row in rows:
        hall_number = int(row[0])
        setting = row[2] or ""
        lg_number = int(row[4])
        hall_symbol = row[6]
        hm_short = row[7]
        hm_full = row[8]
        arith = LG_TO_ARITH[lg_number]
        cent = centering_for(hall_symbol)
        print(
            f"    LayerHallSymbolEntry::new({hall_number}, {lg_number}, "
            f'{arith}, "{setting}", "{hall_symbol}", "{hm_short}", '
            f'"{hm_full}", {cent}),'
        )
    print("];")


if __name__ == "__main__":
    main()
