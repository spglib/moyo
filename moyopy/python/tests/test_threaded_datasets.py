from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import pytest

from moyopy import (
    Cell,
    CollinearMagneticCell,
    MoyoCollinearMagneticDataset,
    MoyoDataset,
    MoyoLayerDataset,
    MoyoNonCollinearMagneticDataset,
    NonCollinearMagneticCell,
)


def test_dataset_constructors_in_thread_pool(
    wurtzite: Cell,
    rutile_type3: CollinearMagneticCell,
    rutile_type3_noncollinear: NonCollinearMagneticCell,
):
    layer = Cell([[1, 0, 0], [0, 1, 0], [0, 0, 5]], [[0, 0, 0]], [1])
    num_workers = 4
    start = Barrier(num_workers, timeout=30)

    def identify(_: int) -> tuple[int, int, int, int]:
        start.wait()
        return (
            MoyoDataset(wurtzite).number,
            MoyoLayerDataset(layer).number,
            MoyoCollinearMagneticDataset(rutile_type3).uni_number,
            MoyoNonCollinearMagneticDataset(rutile_type3_noncollinear).uni_number,
        )

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        numbers = list(executor.map(identify, range(num_workers)))

    assert numbers == [(186, 61, 1158, 496)] * num_workers


def test_dataset_error_in_thread_pool(wurtzite: Cell):
    tilted_layer = Cell([[1, 0, 0], [0, 1, 0], [0.5, 0, 5]], [[0, 0, 0]], [1])

    with ThreadPoolExecutor(max_workers=1) as executor:
        failed_search = executor.submit(MoyoLayerDataset, tilted_layer)
        with pytest.raises(ValueError, match="Aperiodic axis c is not perpendicular"):
            failed_search.result(timeout=30)

        next_search = executor.submit(MoyoDataset, wurtzite)
        assert next_search.result(timeout=30).number == 186
