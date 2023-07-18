# -*- coding: utf-8 -*-
import pytest


@pytest.fixture
def reads():
    """Reads to use in multiple test cases, peak at read 7"""
    return array(
        [
            [1, 5],
            [2, 6],
            [3, 8],
            [4, 10],
            [5, 8],
            [6, 10],
            [7, 1200],
            [8, 14],
            [9, 1],
            [10, 2],
            [11, 5],
            [12, 6],
            [109, 2],
            [208, 4],
        ]
    )


@pytest.fixture
def partially_empty_reads():
    """Example reads with low density/missing reads.  Peak at 7"""
    return array([[1, 1], [2, 1], [5, 2], [7, 1200], [8, 1], [9, 1], [11, 3], [15, 1]])
