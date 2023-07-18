# -*- coding: utf-8 -*-
import sys

import pytest
from mock import patch
from numpy import append, array
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
    @pytest.mark.parameterize("gap,w_sz", [(1, 10), (5, 50), (100, 100)])
    def test_add_padding_new_length(self, gap, w_sz, partially_empty_reads):
        """Test that padding is correctly added to partially empty reads"""
        new_reads = _add_padding(partially_empty_reads, gap, w_sz, random_seed=1)
        assert len(new_reads[0, :]) == (partially_empty_reads[0, -1] + gap + w_sz) - (
            partially_empty_reads[0, 0] - gap - w_sz
        )

    @pytest.mark.parameterize("gap,w_sz", [(1, 10), (5, 50), (100, 100)])
    def test_add_padding_new_length(self, gap, w_sz, partially_empty_reads):
        """Test that padding is correctly added to partially empty reads"""
        new_reads = _add_padding(partially_empty_reads, gap, w_sz, random_seed=1)
        assert len(new_reads[0, :]) == (partially_empty_reads[0, -1] + gap + w_sz) - (
            partially_empty_reads[0, 0] - gap - w_sz
        )


class TestMainAndParseArgsZscore:
    @pytest.fixture
    def regular_argslist(self):
        """A normal list of sys.argv[1:]"""
        return [
            "test_file",
            "--gap",
            "1",
            "--w_sz",
            "3",
            "--save_file",
            False,
        ]

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
        reads = array([[a, b] for a, b in zip(range(1, 1000), range(1001, 2000))])
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


class TestZScores:
    def test_z_scores_regular(self, reads):
        """Z-scores of the reads fixture"""
        print(z_scores(reads, gap=1, w_sz=3))
        assert_array_almost_equal(
            z_scores(reads, gap=1, w_sz=3),
            array(
                [
                    [5, 0],
                    [6, -0.70262826],
                    [7, 202.20038777],
                    [8, -0.56829815],
                    [9, -0.59611538],
                    [10, -0.58959947],
                ]
            ),
        )

    def test_z_scores_outlier(self, reads):
        """An outlier (near the edge where peaks aren't found) doesn't affect score"""
        reads[11] = [12, 1e8]
        print(z_scores(reads, gap=1, w_sz=3))
        assert_array_almost_equal(
            z_scores(reads, gap=1, w_sz=3),
            array(
                [
                    [5, 0],
                    [6, -0.70262826],
                    [7, 202.20038777],
                    [8, -0.56829815],
                    [9, -0.59611538],
                    [10, -0.58959947],
                ]
            ),
        )


class TestValidateGapWindow:
    def test_window_zero(self):
        """Windows can't be zero"""
        with pytest.raises(ValueError) as e_info:
            _validate_gap_window(100, 0)

        assert (
            e_info.value.args[0]
            == "Window size must be larger than 1 to find a z-score"
        )

    def test_window_gap_positive(self):
        """Windows and gaps can be positive"""
        _validate_gap_window(100, 1)

    def test_gap_negative(self):
        """Gaps can't be negative"""
        with pytest.raises(ValueError) as e_info:
            _validate_gap_window(-1, 100)

        assert (
            e_info.value.args[0] == "Gap size must be at least zero to find a z-score"
        )

    def test_gap_zero(self):
        """Gaps can be zero, but should warn"""
        with pytest.warns(UserWarning):
            _validate_gap_window(0, 100)
