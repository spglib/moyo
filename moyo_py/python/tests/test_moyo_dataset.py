import moyo


def test_moyo_dataset(wurtzite: moyo.Structure):
    dataset = moyo.MoyoDataset(wurtzite)
    import IPython; IPython.embed(colors='neutral')
