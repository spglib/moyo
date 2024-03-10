import moyo


def test_moyo_dataset(wurtzite: moyo.Structure):
    dataset = moyo.MoyoDataset(wurtzite)
    assert dataset.number == 186
