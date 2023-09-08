# -*- coding: utf-8 -*-
"""Functions to help convert aligned bowtie sam files into 4 Rend-seq tracks."""

import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from file_funcs import write_wig


def peak_adjust(peaks, pos):
    """Check to see if other end of a read is a peak and adjust density.

    peaks = pandas series of track to check for peak location.
    pos = nt position to check if the other end is a peak for.
    """
    pos_ind = np.argwhere(np.abs(np.asarray(peaks) - pos) < 3)
    return 1 if len(pos_ind) == 0 else 0


def bowtie_to_wig(infile, wig_file_prefix="", peak_file=None):
    """Given a .sam file alignment from bowtie, convert bowtie to wig file.

    -infile: path to bowtie file
    -wig_file_prefix: optional path to location to put wigs.  If not provided
        the folder containing the original bowtie alignment will be used.
    - peaks: optional (default = None). Path to a peaks file, which can be used
        eliminate the shadows from peaks.
    """
    peaks = pd.read_csv(peak_f) if peak_file else None
    if wig_file_prefix == "":
        wig_file_prefix = os.path.basename(infile)

    genomes = {}

    with open(infile, "r") as f:
        for line in f:
            fields = line.split("\t")
            length = len(fields[9])
            if (length > 14) and (length < 45):
                fiveprime = int(fields[3])
                strand = str(fields[1])
                chrom = str(fields[2])
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
        write_wig(genomes[chrom]["3f"], "".join([wig_file_prefix, "_3f.wig"]), chrom)
        write_wig(genomes[chrom]["3r"], "".join([wig_file_prefix, "_3r.wig"]), chrom)
        write_wig(genomes[chrom]["5f"], "".join([wig_file_prefix, "_5f.wig"]), chrom)
        write_wig(genomes[chrom]["5r"], "".join([wig_file_prefix, "_5r.wig"]), chrom)


if __name__ == "__main__":
    infile = sys.argv[1]
    wig_file_prefix = sys.argv[2] if len(sys.argv) > 2 else ""
    peaks = None
    if len(sys.argv) > 3:
        peak_f = sys.argv[3]
    bowtie_to_wig(infile, wig_file_prefix, peaks)
