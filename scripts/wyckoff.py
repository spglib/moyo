from copy import deepcopy

import click


@click.command()
@click.argument("wyckoff-input", type=click.File("r"))
def main(wyckoff_input):
    all_data = []

    lines = wyckoff_input.readlines()
    hall_number = None
    for line in lines:
        line = line.strip('\n')
        if line[0].isdigit() or (line == "end of data"):
            if line[0].isdigit():
                hall_number = int(line.split(':')[0])
        else:
            contents = line.split(':')
            if contents[2] == "":
                continue
            multiplicity = int(contents[2])
            letter = contents[3]
            site_symmetry = contents[4]
            coordinates = contents[5].strip('(').strip(')')
            all_data.append((hall_number, multiplicity, letter, site_symmetry, coordinates))

    print(f"const WYCKOFF_DATABASE: [WyckoffPosition; {len(all_data)}] = [")
    for data in all_data:
        hall_number, multiplicity, letter, site_symmetry, coordinates = data
        print(f'    WyckoffPosition::new({hall_number}, {multiplicity}, \'{letter}\', "{site_symmetry}", "{coordinates}"),')
    print("];")


if __name__ == '__main__':
    main()
