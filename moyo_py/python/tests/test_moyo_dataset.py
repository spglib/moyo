import moyo


def test_moyo_dataset(wurtzite: moyo.Structure):
    dataset = moyo.MoyoDataset(wurtzite)
    assert dataset.number == 186
    assert dataset.hall_number == 480
