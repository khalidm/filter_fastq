'''
Module      : Main
Description : Filter sequences in FASTQ format with a given prefix.
Copyright   : (c) Khalid Mahmood, 27 Sep 2018
License     : BSD-3-Clause
Maintainer  : kmahmood@unimelb.edu.au
Portability : POSIX
'''

#!/usr/bin/python

import sys
import os
import argparse
# from argparse import ArgumentParser
import getopt
import sys
import logging
import pkg_resources
from Bio import SeqIO

import pandas as pd
import numpy as np
#from cyvcf2 import VCF
#import pysam

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
DEFAULT_VERBOSE = False
PROGRAM_NAME = "filter_fastq"

# def exit_with_error(message, exit_status):
#     '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
#     Then exit program with supplied exit status.
#     Arguments:
#         message: an error message as a string.
#         exit_status: a positive integer representing the exit status of the
#         program.
#     '''
#     logging.error(message)
#     print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
#     sys.exit(exit_status)

def get_fq_names(fq_to_parse):
    filename_w_ext = os.path.basename(fq_to_parse)
    filename, file_extension = os.path.splitext(filename_w_ext)
    fq_list = filename.split("_")
    return fq_list

def make_fastq_reads(inputfq, primer1, primer2):
    filtered_seqs = []
    for record in SeqIO.parse(inputfq, "fastq"):
        if not (record.seq.startswith(primer1)) or not (record.seq.startswith(primer2)):
            filtered_seqs.append(record)
    return filtered_seqs

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, dest="input", help="Input file containing prefix sequences", required=True)
    parser.add_argument("-f1", "--fastq1", type=str, dest="fq1", help="Input fastq file (FASTQ) R1", required=True)
    parser.add_argument("-f2", "--fastq2", type=str, dest="fq2", help="Input fastq file (FASTQ) R2", required=True)
    parser.add_argument("-o", "--output", type=str, dest="out", help="Output file (tabular)", required=True)
    parser.add_argument("-v", "--verbosity", action="count", default=0)

    args = parser.parse_args()
    input = args.input
    outputfile_name = args.out
    fq1 = args.fq1
    fq2 = args.fq2

    df = pd.read_csv(input, sep=',', header=0, index_col="offtarget_rank_amplicon_name")
    # Sample38_BSTP-S134_L01_R1_001.fastq
    # Sample38_BSTP-S134_L01_R2_001.fastq
    # Sample38_BSTP-S134.csv
    fq1_list = get_fq_names(fq1)
    fq2_list = get_fq_names(fq2)
    for index, row in df.iterrows():
        new_fq1 = fq1_list[0] + "-" + str(index).replace("_", "-") + "_" + fq1_list[1] + "_" + fq1_list[2] + "_" + fq1_list[3] + "_" + fq1_list[4] + ".fastq"
        new_fq2 = fq2_list[0] + "-" + str(index).replace("_", "-") + "_" + fq2_list[1] + "_" + fq2_list[2] + "_" + fq2_list[3] + "_" + fq2_list[4] + ".fastq"
        print "Writing fastq file: " + new_fq1
        SeqIO.write(make_fastq_reads(fq1, row.forward_primer_sequence, row.reverse_primer_sequence), new_fq1, "fastq")
        print "Writing fastq file: " + new_fq2
        SeqIO.write(make_fastq_reads(fq2, row.forward_primer_sequence, row.reverse_primer_sequence), new_fq2, "fastq")
        # print new_fq1 + "\t" + new_fq2 + "\t" + row.forward_primer_sequence
        # print str(index) + "\t" + ','.join(row)

    # filtered_seqs = []
    # for record in SeqIO.parse(fq1, "fastq"):
    #     filtered_seqs.append(record)
    #
    # SeqIO.write(filtered_seqs, outputfile_name, "fastq")


if __name__ == "__main__":
    main(sys.argv)
