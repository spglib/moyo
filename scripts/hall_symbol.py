import click
import pandas as pd


def number_to_arithmetic_crystal_class():
    # Ordered the same as https://dictionary.iucr.org/Arithmetic_crystal_class
    # fmt: off
    numbers = list(range(1, 230 + 1))
    arithmetic_numbers = [
        ## Triclinic
        # 1P
        1,
        # -1P
        2,
        ## Monoclinic
        # 2P
        3,  3,
        # 2C
        4,
        # mP
        5,  5,
        # mC
        6,  6,
        # 2/mP
        7,  7,
        # 2/mC
        8,
        # 2/mP
        7,  7,
        # 2/mC
        8,
        ## Orthorhombic
        # 222P
        9,  9, 9,  9,
        # 222C
        10, 10,
        # 222F
        11,
        # 222I
        12, 12,
        # mm2P
        13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
        # mm2C
        14, 14, 14,
        # 2mmC
        15, 15, 15, 15,
        # mm2F
        16, 16,
        # mm2I
        17, 17, 17,
        # mmmP
        18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
        # mmmC
        19, 19, 19, 19, 19, 19,
        # mmmF
        20, 20,
        # mmmI
        21, 21, 21, 21,
        ## Tetragonal
        # 4P
        22, 22, 22, 22,
        # 4I
        23, 23,
        # -4P
        24,
        # -4I
        25,
        # 4/mP
        26, 26, 26, 26,
        # 4/mI
        27, 27,
        # 422P
        28, 28, 28, 28, 28, 28, 28, 28,
        # 422I
        29, 29,
        # 4mmP
        30, 30, 30, 30, 30, 30, 30, 30,
        # 4mmI
        31, 31, 31, 31,
        # -42mP
        32, 32, 32, 32,
        # -4m2P
        33, 33, 33, 33,
        # -4m2I
        34, 34,
        # -42mI
        35, 35,
        # 4/mmmP
        36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36,
        # 4/mmmI
        37, 37, 37, 37,
        ## Trigonal
        # 3P
        38, 38, 38,
        # 3R
        39,
        # -3P
        40,
        # -3R
        41,
        # 312P
        42,
        # 321P
        43,
        # 312P
        42,
        # 321P
        43,
        # 312P
        42,
        # 321P
        43,
        # 32R
        44,
        # 3m1P
        45,
        # 31mP
        46,
        # 3m1P
        45,
        # 31mP
        46,
        # 3mR
        47, 47,
        # -31mP
        48, 48,
        # -3m1P
        49, 49,
        # -3mR
        50, 50,
        ## Hexagonal
        # 6P
        51, 51, 51, 51, 51, 51,
        # -6P
        52,
        # 6/mP
        53, 53,
        # 622P
        54, 54, 54, 54, 54, 54,
        # 6mmP
        55, 55, 55, 55,
        # -6m2P
        57, 57,
        # -62mP
        56, 56,
        # 6/mmmP
        58, 58, 58, 58,
        # 23P
        59,
        # 23F
        60,
        # 23I
        61,
        # 23P
        59,
        # 23I
        61,
        # m-3P
        62, 62,
        # m-3F
        63, 63,
        # m-3I
        64,
        # m-3P
        62,
        # m-3I
        64,
        # 432P
        65, 65,
        # 432F
        66, 66,
        # 432I
        67,
        # 432P
        65, 65,
        # 432I
        67,
        # -43mP
        68,
        # -43mF
        69,
        # -43mI
        70,
        # -43mP
        68,
        # -43mF
        69,
        # -43mI
        70,
        # m-3mP
        71, 71, 71, 71,
        # m-3mF
        72, 72, 72, 72,
        # m-3mI
        73, 73,
    ]
    # fmt: on

    df = pd.DataFrame({"number": numbers, "arithmetic_number": arithmetic_numbers})
    return df

@click.command()
@click.argument("spg_input", type=click.File("r"))
def main(spg_input):
    names = [
        "hall_number",
        "setting",
        "number",
        "hall_symbol",
        "HM_symbol_short_all",
        "HM_symbol_full",
    ]
    df_hall = pd.read_csv(
        spg_input,
        header=None,
        names=names,
        usecols=[
            0,  # hall_number
            2,  # setting
            4,  # number
            6,  # hall symbol
            7,  # HM symbol (short)
            8,  # HM symbol (full)
        ],
    )
    df_hall["setting"].fillna(value="", inplace=True)
    df_hall["HM_symbol_short"] = (
        df_hall.HM_symbol_short_all.str.split("=")
        .apply(lambda lst: lst[0])
        .str.rstrip(" ")
    )
    df_hall.drop(columns=["HM_symbol_short_all"], inplace=True)
    df_hall["centering"] = df_hall['hall_symbol'].apply(lambda s: f"Centering::{s[0]}" if s[0] != '-' else f"Centering::{s[1]}")

    df_arith = number_to_arithmetic_crystal_class()

    df = pd.merge(df_hall, df_arith, on="number", how="left")

    for _, row in df.iterrows():
        print(f'HallSymbolEntry::new({row["hall_number"]}, {row["number"]}, {row["arithmetic_number"]}, \"{row["setting"]}\", \"{row["hall_symbol"]}\", \"{row["HM_symbol_short"]}\", \"{row["HM_symbol_full"]}\", {row["centering"]}),')


if __name__ == "__main__":
    main()
