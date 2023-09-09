# -*- coding: utf-8 -*-

import os
import random
import tempfile
from distutils import dir_util

import pytest
from numpy.testing import assert_array_equal

from rendseq.file_funcs import open_wig
from rendseq.make_wig import bowtie_to_wig


@pytest.fixture
def datadir(tmpdir, request):
    """Will find the file in the corresponding data dir and move it to tmpdir"""
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)
    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))
    return tmpdir


def test_wig_trims_untemplated_addition(datadir):
    temp_dir = tempfile.gettempdir()
    rand_str = "".join(
        random.choice("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(10)
    )
    wig_prefix = "".join([str(temp_dir), "/wigs", rand_str])
    bowtie_file = datadir.join("test_files/test.fastq")
    bowtie_to_wig(bowtie_file, wig_file_prefix=wig_prefix)
    wig_ends = ["3f.wig", "3r.wig", "5f.wig", "5r.wig"]
    for w in wig_ends:
        target_wig_name = "".join(["test_files/test", w])
        target_wig_file = datadir.join(target_wig_name)
        target_w, target_c = open_wig(target_wig_file)
        gen_w, gen_c = open_wig("".join([wig_prefix, w]))
        assert gen_c == target_c
        assert_array_equal(target_w, gen_w)
