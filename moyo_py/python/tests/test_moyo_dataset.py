import moyo


def test_moyo_dataset(wurtzite: moyo.Cell):
    dataset = moyo.MoyoDataset(wurtzite)
    assert dataset.number == 186
    assert dataset.hall_number == 480
