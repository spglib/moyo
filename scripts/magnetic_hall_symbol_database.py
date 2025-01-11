from pathlib import Path

import click
from magnetic_space_group_type import get_msg_numbers, get_msg_table


@click.command()
@click.argument("spglib_dir", type=str)
def main(spglib_dir: str):
    msg_dataset_dir = Path(spglib_dir) / "database" / "msg"
    msg_table = get_msg_table(msg_dataset_dir / "magnetic_hall_symbols.yaml")

    msg_numbers = get_msg_numbers(msg_dataset_dir / "msg_numbers.csv")
    bns_to_uni = {bns_number: uni_number for _, bns_number, _, uni_number in msg_numbers}

    contents = []
    for hall_number, dct in msg_table.items():
        for bns_number, mhall_symbol in dct.items():
            uni_number = bns_to_uni[bns_number]
            mhall_symbol = mhall_symbol.replace('"', "=")
            line = f'MagneticHallSymbolEntry::new("{mhall_symbol}", {uni_number}),'
            contents.append(line)

    print("const MAGNETIC_HALL_SYMBOL_DATABASE: [MagneticHallSymbolEntry; 1651] = [")
    for line in contents:
        print(f"    {line}")
    print("];")


if __name__ == "__main__":
    main()
