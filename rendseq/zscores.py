# -*- coding: utf-8 -*-
"""Functions needed for z-score transforming raw rendSeq data."""
import argparse
import sys
import warnings
from os.path import abspath

import numpy as np
from numpy import mean, nan, std, zeros

from rendseq.file_funcs import _validate_reads, make_new_dir, open_wig, write_wig


def _get_lower_reads(cur_ind, target_start, target_stop, reads):
    """Calculate the padded lower reads index in range for the z-score calculation."""
    cur_ind = max(min(cur_ind, len(reads) - 1), 0)
    vals = np.random.normal(0, 1, [abs(target_stop - target_start) + 1, 1])
    ind = 0
    while cur_ind >= 0 and reads[cur_ind, 0] >= target_stop:
        vals[ind] = reads[cur_ind, 1]
        cur_ind -= 1
        ind += 1
    return vals


def _get_upper_reads(cur_ind, target_start, target_stop, reads):
    """Fetch the padded upper reads needed for z score calculation with zero padding."""
    cur_ind = min(max(cur_ind, 0), len(reads) - 1)

    vals = np.random.normal(0, 1, [abs(target_stop - target_start) + 1, 1])
    ind = 0
    while reads[cur_ind, 0] < target_stop and not cur_ind >= len(reads) - 1:
        vals[ind] = reads[cur_ind, 1]
        ind += 1
        cur_ind += 1
    return vals


def _calc_z_score(vals, calc_val):
    """Calculate a z-score given a value, mean, and standard deviation.

    vals = values to calculate the z score with respect to.
    calc_val = value to calculate z score for.

    NOTE: The z_score() of a constant vector is 0
    """
    v_std = std(vals)
    v_mean = mean(vals)
    if nan in [v_mean, v_std]:
        return calc_val
    if v_std == 0:
        return 0 if v_mean == calc_val else (calc_val - v_mean) / 0.2
    return (calc_val - v_mean) / v_std


def _remove_outliers(vals, method="remove_by_std"):
    """Normalize window of reads by removing outliers (values 2.5 std > mean).

    Parameters
    ----------
        -vals: an array of raw read values to be processed

    Returns
    -------
        -new_v: another array of raw values which has had the extreme values
            removed.
    """
    normalized_vals = vals

    if method == "remove_by_std" and len(vals) > 1:
        v_std = std(vals)
        if v_std != 0:
            normalized_vals = [v for v in vals if abs(_calc_z_score(vals, v)) < 2.5]

    return normalized_vals


def _validate_gap_window(gap, w_sz):
    """Check that gap and window size are reasonable in r/l_score_helper."""
    if w_sz < 1:
        raise ValueError("Window size must be larger than 1 to find a z-score")
    if gap < 0:
        raise ValueError("Gap size must be at least zero to find a z-score")
    if gap == 0:
        warnings.warn(
            "Warning...a gap size of 0 includes the current position.", stacklevel=2
        )


def z_scores(reads, gap=5, w_sz=50):
    """Perform modified z-score transformation of reads.

    Parameters
    ----------
        -reads 2xn array - raw rendseq reads
        -gap (interger):   number of reads surround the current read of
            interest that should be excluded in the z_score calculation.
        -w_sz (integer): the max distance (in nt) away from the current position
            one should include in zscore calulcation.

    Returns
    -------
        -z_scores (2xn array): a 2xn array with the first column being position
            and the second column being the z_score.
    """
    _validate_gap_window(gap, w_sz)
    _validate_reads(reads)
    # make array of zscores - same length as raw reads, trimming based on window size:
    z_scores = zeros([len(reads) - 2 * (gap + w_sz), 2])

    # first column of return array is the location of the raw reads
    z_scores[:, 0] = reads[gap + w_sz : len(reads) - (gap + w_sz), 0]

    # Iterate through each valid read, recording z-score
    for i in range((gap + w_sz + 1), (len(reads) - (gap + w_sz))):
        # calculate the z score with values from the left:
        i_score_pos = i - (gap + w_sz)
        if reads[i, 1] > 5:
            l_vals = _get_upper_reads(
                i + gap, reads[i, 0] + gap, reads[i, 0] + gap + w_sz, reads
            )
            l_score = _calc_z_score(_remove_outliers(l_vals), reads[i, 1])

            # calculate z score with reads from the right:
            r_vals = _get_lower_reads(
                i - gap, reads[i, 0] - gap, reads[i, 0] - gap - w_sz, reads
            )
            r_score = _calc_z_score(_remove_outliers(r_vals), reads[i, 1])

            # set the zscore to be the smaller valid score of the left/right scores.
            z_scores[i_score_pos, 1] = (
                r_score if abs(r_score) < abs(l_score) else l_score
            )
        else:
            z_scores[i_score_pos, 1] = 1

    return z_scores


def parse_args_zscores(args):
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
    args = parse_args_zscores(sys.argv[1:])

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
