# -*- coding: utf-8 -*-
import sys

import pytest
from mock import patch
from numpy import array, concatenate, mean, std
from numpy.random import normal
from numpy.testing import assert_array_almost_equal, assert_array_equal

from rendseq.make_peaks import thresh_peaks
from rendseq.zscores import z_scores

step_noise_len = 100
step_internal_len = 100


@pytest.fixture
def step_up_peak():
    """Make a data array that steps up from the noise floor with peak."""
    noise = array(
        [
            [loc, z]
            for loc, z in zip(
                range(1, step_noise_len + 1),
                normal(0, 0.2, size=(1, step_noise_len))[0],
            )
        ]
    )
    peak = [[step_noise_len + 1, 40]]
    internal = array(
        [
            [loc, z]
            for loc, z in zip(
                range(step_noise_len + 2, step_noise_len + step_internal_len),
                normal(5, 2, size=(1, step_internal_len))[0],
            )
        ]
    )
    return concatenate((noise, peak, internal), axis=0)


@pytest.fixture
def step_down_peak():
    """Make a data array that steps up from the noise floor with peak."""
    internal = array(
        [
            [loc, z]
            for loc, z in zip(
                range(1, step_internal_len + 1),
                normal(5, 2, size=(1, step_internal_len))[0],
            )
        ]
    )
    peak = [[step_internal_len + 1, 100]]
    noise = array(
        [
            [loc, z]
            for loc, z in zip(
                range(step_internal_len + 2, step_noise_len + step_internal_len),
                normal(5, 2, size=(1, step_noise_len))[0],
            )
        ]
    )
    return concatenate((internal, peak, noise), axis=0)


def test_thresh_step_up(step_up_peak):
    """Test we can call peaks that step up from the noise floor."""
    z_scr = z_scores(step_up_peak, gap=1, w_sz=5)
    peaks = thresh_peaks(z_scr, thresh=5)
    called_peaks = [p[1] for p in peaks]
    peak_arr = (
        [0 for _ in range(step_noise_len - 6)]
        + [1]
        + [0 for _ in range(step_internal_len - 8)]
    )
    assert_array_equal(called_peaks, peak_arr)


def test_thresh_step_down(step_down_peak):
    """Test we can call peaks that step down to the noise floor."""
    z_scr = z_scores(step_down_peak, gap=1, w_sz=5)
    peaks = thresh_peaks(z_scr, thresh=5)
    called_peaks = [p[1] for p in peaks]
    peak_arr = (
        [0 for _ in range(step_internal_len - 6)]
        + [1]
        + [0 for _ in range(step_noise_len - 8)]
    )
    assert_array_equal(called_peaks, peak_arr)
