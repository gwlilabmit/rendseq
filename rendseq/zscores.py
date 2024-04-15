# -*- coding: utf-8 -*-
"""Functions needed for z-score transforming raw rendSeq data."""
import argparse
import datetime
import sys
import warnings
from os import mkdir
from os.path import abspath, basename, dirname, exists

import numpy as np

from rendseq.file_funcs import _validate_reads, make_new_dir, open_wig, write_wig


def _add_padding(reads, gap, w_sz):
    """Add gaussian padding to the parts of the original array missing values."""
    start = int(reads[0, 0] - gap - w_sz)
    stop = int(reads[-1, 0] + gap + w_sz)
    padded_reads = np.zeros([stop - start, 2])
    padded_reads[:, 0] = list(range(start, stop))
    padded_reads[:, 1] = np.random.normal(0, 1, stop - start)
    padded_reads[(reads[:, 0].astype(int) - int(reads[0, 0]) + gap + w_sz), 1] = reads[
        :, 1
    ]
    return padded_reads


def _apply_std_dis_cutoff(window):
    new_window = window.copy()
    mask = np.abs(window - np.mean(window)) > 1.5 * np.std(window)
    new_window[mask] = np.mean(window[~mask])
    return new_window


def _get_means_sds(reads, w_sz, percent_trim=0.01, winsorize=True):
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
            padded_reads[gap + w_sz : -(gap + w_sz), 1],
            means[(2*gap + w_sz + 1) : ],
        ),
        sds[(2*gap + w_sz + 1) : ],
    )
    lower_zscores = np.divide(
        np.subtract(
            padded_reads[gap + w_sz : - (gap + w_sz), 1],
            means[: - (w_sz + 2*gap + 1)],
        ),
        sds [: - (w_sz + 2*gap + 1)],
    )
    zscores = padded_reads[gap + w_sz : pad_len - (gap + w_sz)].copy()
    zscores[:, 1] = np.min([lower_zscores, upper_zscores], axis=0)

    return zscores[np.where(zscores[:, 1] > 5), :][0]


def call_all_zscores(wig_prefix, z_score_prefix="", wig_ends=None, zscore_ends=None):
    """Call all z-score for all wigs with wig_prefix and wig_ends.

    Parameters
    ----------
        -wig_prefix: (required), str with the start of the file path for
            all wig files.
        -z_score_prefix: (optional, default=""), name of z_score prefix to use for
            all generated z_score files.  If none is provided one will be generated
            in the directory above the wig directory called zscores.
        -wig_ends: (optional, default=None) list of wig ends to append to
            wig_prefix to make all the wig files for processing.
        -zscore_ends: (optional, default=None) list the z_score ends to use for each
            corresponding wig_end.

    Returns
    -------
        -z_score_prefix: str of the z_score prefix of all z_scores generated here.
    """
    if z_score_prefix == "":
        p_dir = "".join([dirname(dirname(wig_prefix)), "/zscores/"])
        file_b = basename(wig_prefix)
        if not exists(p_dir):
            mkdir(p_dir)
        current_date = str(datetime.date.today())
        z_score_prefix = "".join([p_dir, file_b, "_", current_date])

    if wig_ends is None:
        wig_ends = ["_5f.wig", "_3f.wig", "_5r.wig", "_3r.wig"]
        zscore_ends = [
            "_5f_zscore.wig",
            "_3f_zscore.wig",
            "_5r_zscore.wig",
            "_3r_zscore.wig",
        ]
    for ind, w in enumerate(wig_ends):
        reads, chrom = open_wig("".join([wig_prefix, w]))
        z_scr = z_scores(reads, percent_trim=0.02)
        write_wig(z_scr, "".join([z_score_prefix, zscore_ends[ind]]), chrom)
    return z_score_prefix


def call_all_zscores(wig_prefix, z_score_prefix="", wig_ends=None, zscore_ends=None):
    """Call all z-score for all wigs with wig_prefix and wig_ends.

    Parameters
    ----------
        -wig_prefix: (required), str with the start of the file path for
            all wig files.
        -z_score_prefix: (optional, default=""), name of z_score prefix to use for
            all generated z_score files.  If none is provided one will be generated
            in the directory above the wig directory called zscores.
        -wig_ends: (optional, default=None) list of wig ends to append to
            wig_prefix to make all the wig files for processing.
        -zscore_ends: (optional, default=None) list the z_score ends to use for each
            corresponding wig_end.

    Returns
    -------
        -z_score_prefix: str of the z_score prefix of all z_scores generated here.
    """
    if z_score_prefix == "":
        p_dir = "".join([dirname(dirname(wig_prefix)), "/zscores/"])
        file_b = basename(wig_prefix)
        if not exists(p_dir):
            mkdir(p_dir)
        z_score_prefix = "".join([p_dir, file_b])

    if wig_ends is None:
        wig_ends = ["_5f.wig", "_3f.wig", "_5r.wig", "_3r.wig"]
        zscore_ends = [
            "_5f_zscore.wig",
            "_3f_zscore.wig",
            "_5r_zscore.wig",
            "_3r_zscore.wig",
        ]
    for ind, w in enumerate(wig_ends):
        reads, chrom = open_wig("".join([wig_prefix, w]))
        z_scr = z_scores(reads, percent_trim=0.02)
        write_wig(z_scr, "".join([z_score_prefix, zscore_ends[ind]]), chrom)
    return z_score_prefix


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
