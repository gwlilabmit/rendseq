# -*- coding: utf-8 -*-
"""Functions needed for z-score transforming raw rendSeq data."""
import argparse
import sys
import warnings
from os.path import abspath

import numpy as np

from rendseq.file_funcs import _validate_reads, make_new_dir, open_wig, write_wig


def _add_padding(reads, gap, w_sz):
    """Add gaussian padding to the parts of the original array missing values."""
    start = reads[0, 0] - gap - w_sz
    stop = reads[-1, 0] + gap + w_sz
    padded_reads = np.zeros([stop - start, 2])
    padded_reads[:, 0] = list(range(start, stop))
    padded_reads[:, 1] = np.random.normal(0, 1, stop - start)
    padded_reads[(reads[:, 0] - reads[0, 0] + gap + w_sz), 1] = reads[:, 1]
    return padded_reads


def _apply_std_dis_cutoff(window):
    mask = np.abs(window - np.mean(window)) > 2 * np.std(window)
    return window[~mask]


def _get_means_sds(reads, w_sz, percent_trim=0, winsorize=False):
    """Calculate the arrays of mean and sd of tiled windows along data."""
    if percent_trim > 0:
        sliding_windows = np.sort(
            np.lib.stride_tricks.sliding_window_view(reads[:, 1], w_sz), axis=1
        )[
            :, : int(w_sz * (1 - percent_trim))
        ]  # Note - remove the top ~x% of values (other peaks?)
    else:
        sliding_windows = np.lib.stride_tricks.sliding_window_view(reads[:, 1], w_sz)
    if winsorize:
        sliding_windows = np.apply_along_axis(_apply_std_dis_cutoff, 1, sliding_windows)
    means = np.mean(sliding_windows, axis=1)
    sds = np.std(sliding_windows, axis=1)

    return means, sds


def _validate_gap_window(gap, w_sz):
    """Check that gap and window size are reasonable in r/l_score_helper."""
    if int(w_sz * 0.8) < 1:
        raise ValueError("Window size must be larger than 1 to find a z-score")
    if gap < 0:
        raise ValueError("Gap size must be at least zero to find a z-score")
    if gap == 0:
        warnings.warn(
            "Warning...a gap size of 0 includes the current position.", stacklevel=2
        )


def z_scores(reads, gap=5, w_sz=50, percent_trim=0, winsorize=True):
    """Perform modified z-score transformation of reads.

    Parameters
    ----------
        -reads 2xn array - raw rendseq reads
        -gap (interger):   number of reads surround the current read of
            interest that should be excluded in the z_score calculation.
        -w_sz (integer): the max distance (in nt) away from the current position
            one should include in zscore calulcation.
        -percent_trim - what fraction of the top reads should be dropped before
            calculating the mean and std?  ie 0.1 means the top 10% of reads
            are dropped.
        -winsorize - bool for whether or not after trimming any reads more than
            1.5 std from the mean should be dropped.

    Returns
    -------
        -z_scores (2xn array): a 2xn array with the first column being position
            and the second column being the z_score.
    """
    _validate_gap_window(gap, w_sz)
    _validate_reads(reads)
    padded_reads = _add_padding(reads, gap, w_sz)
    pad_len = len(padded_reads[:, 0])
    means, sds = _get_means_sds(padded_reads, w_sz, percent_trim, winsorize)
    upper_zscores = np.divide(
        np.subtract(
            padded_reads[gap + w_sz : pad_len - (gap + w_sz), 1],
            means[(gap + w_sz) : len(means) - 1 - gap],
        ),
        sds[(gap + w_sz) : len(means) - 1 - gap],
    )
    lower_zscores = np.divide(
        np.subtract(
            padded_reads[gap + w_sz : pad_len - (gap + w_sz), 1],
            means[gap : len(means) - 1 - (gap + w_sz)],
        ),
        sds[gap : len(means) - 1 - (gap + w_sz)],
    )
    zscores = padded_reads[gap + w_sz : pad_len - (gap + w_sz)].copy()
    zscores[:, 1] = np.min([lower_zscores, upper_zscores], axis=0)

    return zscores


def _parse_args_zscores(args):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Takes raw read file and\
                                        makes a modified z-score for each\
                                        position. Takes several optional\
                                        arguments"
    )
    parser.add_argument(
        "filename",
        help="Location of the raw_reads file that\
                                        will be processed using this function.\
                                        Should be a properly formatted wig\
                                        file.",
    )
    parser.add_argument(
        "--gap",
        help="gap (interger):   number of reads\
                                        surround the current read of interest\
                                        that should be excluded in the z_score\
                                        calculation. Defaults to 5.",
        default=5,
    )
    parser.add_argument(
        "--w_sz",
        help="w_sz (integer): the max dis (in nt)\
                                        away from the current position one\
                                        should include in zscore calulcation.\
                                        Default to 50.",
        default=50,
    )
    parser.add_argument(
        "--save_file",
        help="Save the z_scores file as a new\
                                        wig file in addition to returning the\
                                        z_scores.  Default = True",
        default=True,
    )
    return parser.parse_args(args)


def main_zscores():
    """Run Z-score calculations.

    Effect: Writes messages to standard out. If --save-file flag,
    also writes output to disk.
    """
    args = _parse_args_zscores(sys.argv[1:])

    # Calculate z-scores
    filename = args.filename
    print(f"Calculating zscores for file {filename}.")
    reads, chrom = open_wig(filename)
    z_score = z_scores(reads, gap=int(args.gap), w_sz=int(args.w_sz))

    # Save file, if applicable
    if args.save_file:
        _save_zscore(filename, z_score, chrom)
    print(
        "\n".join(
            [
                "Ran zscores.py with the following settings:",
                f"gap: {args.gap}, w_sz: {args.w_sz},",
                f"file_name: {args.filename}",
            ]
        )
    )


def _save_zscore(filename, z_score, chrom):
    filename = abspath(filename).replace("\\", "/")
    file_loc = filename[: filename.rfind("/")]
    z_score_dir = make_new_dir([file_loc, "/Z_scores"])
    file_start = filename[filename.rfind("/") : filename.rfind(".")]
    z_score_file = "".join([z_score_dir, file_start, "_zscores.wig"])
    write_wig(z_score, z_score_file, chrom)
    print(f"Wrote z_scores to {z_score_file}")


if __name__ == "__main__":
    main_zscores()
