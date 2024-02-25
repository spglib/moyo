import pytest
import moyo


def test_sum_as_string():
    assert moyo.sum_as_string(1, 1) == "2"
