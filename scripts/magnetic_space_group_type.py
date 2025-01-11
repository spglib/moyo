import csv
import re
from pathlib import Path
from typing import Any

import click
from ruamel.yaml import YAML

MAX_BNS_NUMBER_LENGTH = 7
MAX_OG_NUMBER_LENGTH = 11


# hall_number -> bns_number -> magnetic Hall symbol
def get_msg_table(path: Path) -> dict[str, dict[str, str]]:
    # Load MSG for ITA standard settings
    with open(path) as f:
        all_datum = dict(YAML().load(f))
    return all_datum


def get_msg_numbers(path: Path) -> list[dict[str, Any]]:
    all_datum = []
    with open(path) as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)  # skip header
        for row in reader:
            if len(row) == 0:
                break

            litvin_number, bns_number, og_number, uni_number = row
            all_datum.append(
                (
                    int(litvin_number),
                    bns_number,
                    og_number,
                    int(uni_number),
                )
            )

    assert len(all_datum) == 1651
    return all_datum


def get_type_of_msg(hall_symbol: str) -> int:
    # awkward classification...
    if "'" not in hall_symbol:
        return 1  # type-I
    if " 1'" in hall_symbol:
        return 2  # type-II
    if len(re.findall(r" 1[a-z]+'", hall_symbol)) > 0:
        return 4  # type-IV
    return 3  # type-III


@click.command()
@click.argument("spglib_dir", type=str)
def main(spglib_dir: str):
    msg_dataset_dir = Path(spglib_dir) / "database" / "msg"
    msg_dataset = get_msg_numbers(msg_dataset_dir / "msg_numbers.csv")

    msg_table = get_msg_table(msg_dataset_dir / "magnetic_hall_symbols.yaml")

    bns_to_mhall = {}
    for dct in msg_table.values():
        bns_to_mhall.update(dct)

    contents = []
    construct_type_counts = {t: 0 for t in [1, 2, 3, 4]}
    for i, (litvin_number, bns_number, og_number, uni_number) in enumerate(msg_dataset):
        assert uni_number == i + 1
        assert 0 < len(bns_number) <= MAX_BNS_NUMBER_LENGTH
        assert 0 < len(og_number) <= MAX_OG_NUMBER_LENGTH

        magnetic_hall_symbol = bns_to_mhall[bns_number]
        construct_type = get_type_of_msg(magnetic_hall_symbol)
        xsg_number = int(bns_number.split(".")[0])
        construct_type_counts[construct_type] += 1

        line = f'MagneticSpaceGroupType::new({uni_number}, {litvin_number}, "{bns_number}", "{og_number}", {xsg_number}, ConstructType::Type{construct_type}),'  # noqa: E501
        contents.append(line)

    # Sanity check
    assert construct_type_counts[1] == 230
    assert construct_type_counts[2] == 230
    assert construct_type_counts[3] == 674
    assert construct_type_counts[4] == 517

    print("const MAGNETIC_SPACE_GROUP_TYPES: [MagneticSpaceGroupType; 1651] = [")
    for line in contents:
        print(f"    {line}")
    print("];")


if __name__ == "__main__":
    main()
