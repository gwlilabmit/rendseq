# -*- coding: utf-8 -*-
"""Functions to help convert aligned bowtie sam files into 4 Rend-seq tracks."""

import json
import os
import sys
from collections import defaultdict

import numpy as np

from rendseq.file_funcs import write_wig


def peak_adjust(peaks, pos):
    """Check to see if other end of a read is a peak and adjust density.

    Parameters
    ----------
        - peaks = pandas series of track to check for peak location.
        - pos = nt position to check if the other end is a peak for.
    """
    pos_ind = np.argwhere(np.abs(np.asarray(peaks) - pos) < 3)
    return 1 if len(pos_ind) == 0 else 0


def write_wrap(genomes, chrom, track, prefix):
    """Write provided data to the wig file.

    Parameters
    ----------
        - genomes: a nested dict structure with 4 tracks per chrom.
        - chrom: the chromosome to write.
        - track: the track to write (3f, 5r etc)
        - prefix: the prefix of the wig file.
    """
    wig_file = "".join([prefix, track, ".wig"])
    data_track = genomes[chrom][track]
    wig_array = np.asarray(
        list(map(lambda x: [x, data_track[x]], sorted(data_track.keys())))
    )
    write_wig(wig_array, wig_file, chrom)


def bowtie_to_wig(infile, wig_file_prefix="", peak_file=None):
    """Given a .sam file alignment from bowtie, convert bowtie to wig file.

    Parameters
    ----------
        - infile: path to bowtie file
        -wig_file_prefix: optional path to location to put wigs.  If not provided
            the folder containing the original bowtie alignment will be used.
        - peaks: optional (default = None). Path to a peaks file, which can be used
            eliminate the shadows from peaks.
    """
    peaks = None
    if peak_file:
        with open(peak_file, "r") as f:
            peaks = json.load(f)

    if wig_file_prefix == "":
        dir = os.path.dirname(infile)
        file = os.path.basename(infile)
        if peaks:
            wig_file_prefix = "".join(
                [dir, "/wigs/", file[: file.rfind(".")]], "_no_shadow"
            )
        else:
            wig_file_prefix = "".join([dir, "/wigs/", file[: file.rfind(".")]])
    genomes = {}

    with open(infile, "r") as f:
        for line in f:
            fields = line.split("\t")
            length = len(fields[9])
            if (length > 14) and (length < 45):
                fiveprime, strand, chrom = (
                    int(fields[3]),
                    str(fields[1]),
                    str(fields[2]),
                )
                strand = "+" if strand == "0" else "-"
                mismatch_info = fields[12][fields[12].rfind(":") + 1 :]
                if strand == "+" and mismatch_info[0] == "0":
                    fiveprime += 1
                    length -= 1
                if strand == "-" and mismatch_info[-1] == "0":
                    length -= 1

                end_5 = fiveprime + 1
                end_3 = fiveprime + length
                if chrom not in genomes.keys():
                    genomes[chrom] = defaultdict(lambda: defaultdict(int))
                if strand == "+":
                    genomes[chrom]["3f"][end_3] += (
                        peak_adjust(peaks["5f"], end_5) if peaks else 1
                    )
                    genomes[chrom]["5f"][end_5] += (
                        peak_adjust(peaks["3f"], end_3) if peaks else 1
                    )
                elif strand == "-":
                    genomes[chrom]["3r"][end_5] += (
                        peak_adjust(peaks["5r"], end_3) if peaks else 1
                    )
                    genomes[chrom]["5r"][end_3] += (
                        peak_adjust(peaks["3r"], end_5) if peaks else 1
                    )

    for chrom in genomes:
        write_wrap(genomes, chrom, "3f", wig_file_prefix)
        write_wrap(genomes, chrom, "3r", wig_file_prefix)
        write_wrap(genomes, chrom, "5f", wig_file_prefix)
        write_wrap(genomes, chrom, "5r", wig_file_prefix)

    return wig_file_prefix


if __name__ == "__main__":
    infile = sys.argv[1]
    wig_file_prefix = sys.argv[2] if len(sys.argv) > 2 else ""
    peaks = sys.argv[3] if len(sys.argv) > 3 else None
    bowtie_to_wig(infile, wig_file_prefix, peaks)
