# -*- coding: utf-8 -*-
"""Functions needed for z-score transforming raw rendSeq data."""
import warnings

from numpy import mean, std, zeros

from rendseq.file_funcs import validate_reads


def _adjust_down(cur_ind, target_val, reads):
    """Calculate the lower reads index in range for the z-score calculation."""
    validate_reads(reads)
    cur_ind = min(cur_ind, len(reads) - 1)
    while reads[cur_ind, 0] > target_val:
        cur_ind -= 1

        if cur_ind == 0:
            break
    return cur_ind


def _adjust_up(cur_ind, target_val, reads):
    """Calculate the higher reads index in range for the z-score calculation."""
    if len(reads) < 1:
        raise ValueError("requires non-empty reads")

    cur_ind = min(max(cur_ind, 0), len(reads))

    while reads[cur_ind, 0] < target_val:
        if cur_ind >= len(reads) - 1:
            break
        cur_ind += 1

    return cur_ind


def _z_score(val, v_mean, v_std):
    """Calculate a z-score given a value, mean, and standard deviation.

    NOTE: The z_score() of a constant vector is 0
    """
    score = 0 if v_std == 0 else (val - v_mean) / v_std
    return score


def _remove_outliers(vals):
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
    if len(vals) > 1:
        v_mean = mean(vals)
        v_std = std(vals)
        if v_std != 0:
            normalized_vals = [v for v in vals if abs(_z_score(v, v_mean, v_std)) < 2.5]

    return normalized_vals


def _calc_score(vals, min_r, cur_val):
    """Compute the z score.

    Parameters
    ----------
        -vals raw read count values array
        -min_r: the minumum number of reads needed to calculate score
        -cur_val: the value for which the z score is being calculated

    Returns
    -------
        -score: the zscore for the current value, or None if insufficent reads
    """
    score = None
    if sum(vals + cur_val) > min_r:
        v_mean = mean(vals)
        v_std = std(vals)

        score = _z_score(cur_val, v_mean, v_std)

    return score


def score_helper(start, stop, min_r, reads, i):
    """Find the z-score of reads[i] relative to the subsection of reads.

    Goes from start to stop, with a read cutoff of min_r
    """
    reads_outlierless = _remove_outliers(list(reads[start:stop, 1]))
    return _calc_score(reads_outlierless, min_r, reads[i, 1])


def validate_gap_window(gap, w_sz):
    """Check that gap and window size are reasonable in r/l_score_helper."""
    if w_sz < 1:
        raise ValueError("Window size must be larger than 1 to find a z-score")
    if gap < 0:
        raise ValueError("Gap size must be at least zero to find a z-score")
    if gap == 0:
        warnings.warn("Warning...a gap size of 0 includes the current position.")


def _l_score_helper(gap, w_sz, min_r, reads, i):
    """Find the z_score based on reads to the left of the current pos."""
    validate_gap_window(gap, w_sz)
    l_start = _adjust_up(i - (gap + w_sz), reads[i, 0] - (gap + w_sz), reads)
    l_stop = _adjust_up(i - gap, reads[i, 0] - gap, reads)
    return score_helper(l_start, l_stop, min_r, reads, i)


def _r_score_helper(gap, w_sz, min_r, reads, i):
    """Find the z_score based on reads to the right of the current pos."""
    validate_gap_window(gap, w_sz)
    r_start = _adjust_down(i + gap, reads[i, 0] + gap, reads)
    r_stop = _adjust_down(i + gap + w_sz, reads[i, 0] + gap + w_sz, reads)
    return score_helper(r_start, r_stop, min_r, reads, i)


def z_scores(reads, gap=5, w_sz=50, min_r=20):
    """Perform modified z-score transformation of reads.

    Parameters
    ----------
        -reads 2xn array - raw rendseq reads
        -gap (interger):   number of reads surround the current read of
            interest that should be excluded in the z_score calculation.
        -w_sz (integer): the max distance (in nt) away from the current position
            one should include in zscore calulcation.
        -min_r (integer): density threshold. If there are less than this number
            of reads going into the z_score calculation for a point that point
            is excluded.  note this is sum of reads in the window
        -file_name (string): the base file_name, can be passed in to customize
            the message printed

    Returns
    -------
        -z_score (2xn array): a 2xn array with the first column being position
            and the second column being the z_score.
    """
    # make array of zscores - same length as raw reads, trimming based on window size:
    z_score = zeros([len(reads) - 2 * (gap + w_sz), 2])

    # first column of return array is the location of the raw reads
    z_score[:, 0] = reads[gap + w_sz : len(reads) - (gap + w_sz), 0]

    # Iterate through each valid read, recording z-score
    for i in range((gap + w_sz + 1), (len(reads) - (gap + w_sz))):
        # calculate the z score with values from the left:
        l_score = _l_score_helper(gap, w_sz, min_r, reads, i)
        # calculate z score with reads from the right:
        r_score = _r_score_helper(gap, w_sz, min_r, reads, i)

        # The location in which this z-score should go into the final array
        i_score_pos = i - (gap + w_sz)

        # set the zscore to be the smaller valid score of the left/right scores
        # If neither score is valid, Z-score is 0
        z_score[i_score_pos, 1] = 0
        if l_score is not None:
            if r_score is not None:
                z_score[i_score_pos, 1] = (
                    r_score if abs(r_score) < abs(l_score) else l_score
                )
            else:
                z_score[i_score_pos, 1] = l_score

        elif r_score is not None:
            z_score[i_score_pos, 1] = r_score

    return z_score
