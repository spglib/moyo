import moyopy


def test_moyo_dataset(wurtzite: moyopy.Cell):
    dataset = moyopy.MoyoDataset(wurtzite)
    assert dataset.number == 186
    assert dataset.hall_number == 480
