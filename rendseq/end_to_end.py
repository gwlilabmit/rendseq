# -*- coding: utf-8 -*-
"""Functions which tie together rendseq functionality into single call."""
from rendseq.make_peaks import thresh_all_peaks
from rendseq.make_wig import bowtie_to_wig
from rendseq.zscores import call_all_zscores


def bowtie_to_peak(bowtie_file, thresh=12, peak_file=None):
    """Convert bowtie aligned file to peaks for 5r, 5f, 3r, and 3f.

    Parameters
    ----------
        - bowtie_file: (required) str path to the bowtie file.
        - thresh: (optional, default=12) float of what z score threshold
            should be used for thresholding z score values.
        - peaks: (optional, default = None), str path to a peaks file. If
            provided peaks will be used to remove shadows from the peaks.

    Returns
    -------
        -peaks: str path to the location of the peaks file.
    """
    wig_prefix = bowtie_to_wig(bowtie_file, peak_file=peak_file)
    z_score_prefix = call_all_zscores(wig_prefix)
    return thresh_all_peaks(z_score_prefix, thresh=thresh)


def bowtie_to_peak_remove_shadows(bowtie_file, thresh=12):
    """Convert bowtie aligned file to peaks after shadow removal.

    Parameters
    ----------
        - bowtie_file: (required) str path to the bowtie file.

    Returns
    -------
        -peaks: str path to the location of the peaks file (called
            after shadow removal).
    """
    first_pass_peaks = bowtie_to_peak(bowtie_file, thresh=thresh)
    return bowtie_to_peak(bowtie_file, peak_file=first_pass_peaks, thresh=thresh)
