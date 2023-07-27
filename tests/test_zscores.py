# -*- coding: utf-8 -*-
import math
import sys

import numpy as np
import pytest
import scipy.stats as stats
from mock import patch
from numpy.testing import assert_array_almost_equal, assert_array_equal

from rendseq.file_funcs import write_wig
from rendseq.zscores import (
    _add_padding,
    _get_means_sds,
    _parse_args_zscores,
    _save_zscore,
    _validate_gap_window,
    main_zscores,
    z_scores,
)


class TestAddPadding:
    @pytest.mark.parametrize("gap,w_sz", [(1, 10), (5, 50), (100, 100)])
    def test_add_padding_new_length(self, gap, w_sz, partially_empty_reads):
        """Test that padding is correctly added to partially empty reads."""
        new_reads = _add_padding(partially_empty_reads, gap, w_sz)
        assert len(new_reads[:, 0]) == (partially_empty_reads[-1, 0] + gap + w_sz) - (
            partially_empty_reads[0, 0] - gap - w_sz
        )

    @pytest.mark.parametrize("gap,w_sz", [(1, 10), (5, 50), (100, 100)])
    def test_add_padding_is_normal(self, gap, w_sz, partially_empty_reads):
        """Test that added padding is approximately normal."""
        new_reads = _add_padding(partially_empty_reads, gap, w_sz)
        added_indices = list(
            filter(
                lambda x: new_reads[x, 0] not in partially_empty_reads[:, 0],
                list(range(len(new_reads[:, 0]))),
            ),
        )
        added_vals = new_reads[added_indices, 1]
        std = np.std(added_vals)
        mean = np.mean(added_vals)
        sem = std / (len(added_vals) ** (1 / 2))
        assert stats.t.cdf(mean / sem, len(added_vals) - 1) > 0.01


class TestMeanSds:
    @pytest.mark.parametrize("gap,w_sz", [(1, 10), (5, 50), (100, 100)])
    def test_with_no_modification(self, gap, w_sz, reads):
        # Test that with no trimming and no winsorization the calc mean/std is of the unmodified array.
        padded_reads = _add_padding(reads, gap, w_sz)
        calc_mean, calc_std = _get_means_sds(
            padded_reads,
            w_sz,
            percent_trim=0,
            winsorize=False,
        )
        assert calc_mean[0] == np.mean(padded_reads[:w_sz, 1])
        assert calc_std[0] == np.std(padded_reads[:w_sz, 1])

    @pytest.mark.parametrize("gap,w_sz", [(1, 10), (5, 50), (100, 100)])
    def test_trimming(self, gap, w_sz, partially_empty_reads):
        # If we fill the start of the padded array with reads we expect to be trimmed - verify that the
        # mean is just the mean of the array without them.
        padded_reads = _add_padding(partially_empty_reads, gap, w_sz)
        max_in_padded_reads = max(padded_reads[:w_sz, 1])
        padded_reads[: int(w_sz * 0.1)] = max_in_padded_reads * 5
        calc_mean, calc_std = _get_means_sds(
            padded_reads,
            w_sz,
            percent_trim=0.1,
            winsorize=False,
        )
        assert math.isclose(
            calc_mean[0],
            np.mean(padded_reads[int(w_sz * 0.1) : w_sz, 1]),
            rel_tol=1e-07,
        )
        assert math.isclose(
            calc_std[0], np.std(padded_reads[int(w_sz * 0.1) : w_sz, 1]), rel_tol=1e-07
        )

    def test_winsorization(self, winsorization_test_reads):
        # If we fill the start of the padded array with reads we expect to be trimmed - verify that the
        # mean is just the mean of the array without them.
        fake_window = len(winsorization_test_reads[:, 0])
        calc_mean, calc_std = _get_means_sds(
            winsorization_test_reads,
            fake_window,
            percent_trim=0,
            winsorize=True,
        )
        assert math.isclose(
            calc_mean[0],
            np.mean(winsorization_test_reads[1:fake_window, 1]),
            rel_tol=1e-07,
        )
        assert math.isclose(
            calc_std[0],
            np.std(winsorization_test_reads[1:fake_window, 1]),
            rel_tol=1e-07,
        )


class TestMainAndParseArgsZscore:
    def test_main(self, tmpdir, capfd, reads, regular_argslist):
        """Main with normal settings"""
        # Create a wig file with the reads() fixture
        file = tmpdir.join("file.txt")
        chrom = "test_chrom"
        write_wig(reads, file.strpath, chrom)

        # Modify the argslist with the temporary wig file
        regular_argslist = [""] + regular_argslist
        regular_argslist[1] = file.strpath

        # Run main with regular arguments
        with patch.object(sys, "argv", regular_argslist):
            main_zscores()
            out, _ = capfd.readouterr()

        # Expect output
        assert out == "\n".join(
            [
                f"Calculating zscores for file {file.strpath}.",
                f"Ran zscores.py with the following settings:\ngap: 1, w_sz: 3,\nfile_name: {file.strpath}\n",
            ]
        )

    def test_main_defaults(self, tmpdir, capfd):
        """Main with defaults"""
        # Create a wig file with the reads() fixture
        file = tmpdir.join("file.txt")
        chrom = "test_chrom"

        # Need a larger reads file for defaults
        reads = np.array([[a, b] for a, b in zip(range(1, 1000), range(1001, 2000))])
        write_wig(reads, file.strpath, chrom)

        # Modify the argslist with the temporary wig file
        argslist = ["", file.strpath]

        # Run main
        with patch.object(sys, "argv", argslist):
            main_zscores()
            out, err = capfd.readouterr()

        # Expect output
        file_head = file.strpath[:-8].replace("\\", "/")
        assert out == "\n".join(
            [
                f"Calculating zscores for file {file.strpath}.",
                f"Wrote z_scores to {file_head}Z_scores/file_zscores.wig",
                f"Ran zscores.py with the following settings:\ngap: 5, w_sz: 50,\nfile_name: {file.strpath}\n",
            ]
        )

    def test_parse_args(self, regular_argslist):
        """Regular arguments"""
        args = _parse_args_zscores(regular_argslist)
        assert args.filename == "test_file"
        assert args.gap == "1"
        assert args.w_sz == "3"
        assert not args.save_file

    def test_parse_args_defaults(self):
        """Makes sure the arg defaults are as-expected"""
        arg_list = ["test_file"]
        args = _parse_args_zscores(arg_list)

        assert args.filename == "test_file"
        assert args.gap == 5
        assert args.w_sz == 50
        assert args.save_file
