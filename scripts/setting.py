import click
import pandas as pd


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

    print("spglib")
    spglib = []
    for number in range(1, 231):
        df_tmp = df_hall[df_hall["number"] == number].sort_values("hall_number")
        spglib.append(df_tmp.iloc[0]["hall_number"])
    print(spglib)

    print("standard")
    standard = []
    for number in range(1, 231):
        df_tmp = df_hall[df_hall["number"] == number].sort_values("hall_number")
        if df_tmp[df_tmp["setting"] == "2"].empty:
            standard.append(df_tmp.iloc[0]["hall_number"])
        else:
            standard.append(df_tmp[df_tmp["setting"] == "2"].iloc[0]["hall_number"])
    print(standard)


if __name__ == "__main__":
    main()
