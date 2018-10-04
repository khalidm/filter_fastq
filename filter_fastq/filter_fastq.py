'''
Module      : Main
Description : Filter sequences in FASTQ format with a given prefix.
Copyright   : (c) Khalid Mahmood, 27 Sep 2018
License     : BSD-3-Clause
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
# sys.stdout.flush()

from time import gmtime, strftime

import pandas as pd
import numpy as np

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
DEFAULT_VERBOSE = False
PROGRAM_NAME = "filter_fastq"

def get_fq_names(fq_to_parse):
    filename_w_ext = os.path.basename(fq_to_parse)
    filename, file_extension = os.path.splitext(filename_w_ext)
    fq_list = filename.split("_")
    return fq_list

def get_new_fq_names(fq_to_parse):
    filename_w_ext = os.path.basename(fq_to_parse)
    filename, file_extension = os.path.splitext(filename_w_ext)
    fq_list = filename.split("_")
    return fq_list

def make_fastq_reads(inputfq, primer1, primer2):
    filtered_seqs = []
    for record in SeqIO.parse(inputfq, "fastq"):
        if not record.seq.startswith(primer1) and not record.seq.startswith(primer2):
            filtered_seqs.append(record)
    return filtered_seqs

def make_fastq_reads_from_lists(inputfq, primer1, primer2):
    # print "looking for " + str(primer1)
    filtered_seqs = []
    for record in SeqIO.parse(inputfq, "fastq"):
        # print "looking for " + record.seq + " in " + str(primer1) + "/" + str(primer2)
        if not record.seq.startswith(tuple(primer1)) and not record.seq.startswith(tuple(primer2)):
            filtered_seqs.append(record)
    return filtered_seqs

def get_fastq_readids(inputfq, primer1, primer2):
    # print "looking for " + str(primer1)
    filtered_rids = []
    input_seq_iterator = SeqIO.parse(inputfq, "fastq")
    for record in input_seq_iterator:
        #if not record.seq.startswith(tuple(primer1)) and not record.seq.startswith(tuple(primer2)):
        if record.seq.startswith(tuple(primer1)) or record.seq.startswith(tuple(primer2)):
            if record.id not in filtered_rids:
                filtered_rids.append(record.id)
    return filtered_rids

def write_filtered_fastq(filter_list, inputfq, new_fq):
    count = 0
    filtered_seqs = []
    input_seq_iterator = SeqIO.parse(inputfq, "fastq")
    for record in input_seq_iterator:
        if record.id not in filter_list:
        # if any(record.id is not x for x in filter_list):
            # filtered_seqs.append(record)
            yield record
            # c = SeqIO.write(record, new_fq, "fastq")
            # count = count + c
    # return filtered_seqs

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, dest="input", help="Input file containing prefix sequences", required=True)
    parser.add_argument("-f1", "--fastq1", type=str, dest="fq1", help="Input fastq file (FASTQ) R1", required=True)
    parser.add_argument("-f2", "--fastq2", type=str, dest="fq2", help="Input fastq file (FASTQ) R2", required=True)
    # parser.add_argument("-o", "--output", type=str, dest="out", help="Output file (tabular)", required=True)
    parser.add_argument("-v", "--verbosity", action="count", default=0)

    args = parser.parse_args()
    input = args.input
    # outputfile_name = args.out
    fq1 = args.fq1
    fq2 = args.fq2

    df = pd.read_csv(input, sep=',', header=0, index_col="offtarget_rank_amplicon_name")
    forwards = df.forward_primer_sequence.values
    reverse = df.reverse_primer_sequence.values
    # Sample38_BSTP-S134_L01_R1_001.fastq
    # Sample38_BSTP-S134_L01_R2_001.fastq
    # Extract Sample38_BSTP-S134.csv from Sample38_BSTP-S134_L01_R1_001.fastq
    fq1_list = get_fq_names(fq1)
    fq2_list = get_fq_names(fq2)

    # Incremently remove reads
    if len(forwards) == len(reverse):
        for i in range(len(forwards)):
            print "H0 " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
            # print forwards[i] + ":" + reverse[i]
            k = i + 1
            new_fq1 = fq1_list[0] + "-" + str(k) + "_" + fq1_list[1] + "_" + fq1_list[2] + "_" + fq1_list[3] + "_" + fq1_list[4] + ".fastq"
            new_fq2 = fq2_list[0] + "-" + str(k) + "_" + fq2_list[1] + "_" + fq2_list[2] + "_" + fq2_list[3] + "_" + fq2_list[4] + ".fastq"

            filter_ids = {}
            filter_ids_R1 = []
            # filter_ids_R1.append(get_fastq_readids(fq1, forwards[:k], reverse[:k]))
            filter_ids_R1 = get_fastq_readids(fq1, forwards[:k], reverse[:k])
            print "H1 " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
            filter_ids_R2 = []
            # filter_ids_R2.append(get_fastq_readids(fq2, forwards[:k], reverse[:k]))
            filter_ids_R2 = get_fastq_readids(fq2, forwards[:k], reverse[:k])
            print "H2 " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
            # filter_ids = list(set(filter_ids_R1).union(set(filter_ids_R2)))
            filter_ids = set(filter_ids_R1).union(set(filter_ids_R2))
            print "H3 " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

            # input_seq_iterator = SeqIO.parse(, "fastq")
            # short_seq_iterator = (record for record in input_seq_iterator if len(record.seq) < 300)

            # list contianing read_ids to exclude
            # filter_ids
            # seqs = write_filtered_fastq(filter_ids, fq1, new_fq1)
            # print "H3A" + strftime("%Y-%m-%d %H:%M:%S", gmtime())
            # print filter_ids
            c1 = SeqIO.write(write_filtered_fastq(filter_ids, fq1, new_fq1), new_fq1, "fastq")
            print "Writing fastq file: " + new_fq1 + " (" + str(c1) + ")"
            print "H4 " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
            c2 = SeqIO.write(write_filtered_fastq(filter_ids, fq2, new_fq2), new_fq2, "fastq")
            print "Writing fastq file: " + new_fq2 + " (" + str(c2) + ")"
            print "H5 " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

            fq1 = new_fq1
            fq2 = new_fq2
            del filter_ids_R1[:]
            del filter_ids_R2[:]
            filter_ids.clear()


    else:
        print "error message:\n" + parser.print_usage()

    # Iteratively remove only corresponding reads
    # for index, row in df.iterrows():
    #     new_fq1 = fq1_list[0] + "-" + str(index).replace("_", "-") + "_" + fq1_list[1] + "_" + fq1_list[2] + "_" + fq1_list[3] + "_" + fq1_list[4] + ".fastq"
    #     new_fq2 = fq2_list[0] + "-" + str(index).replace("_", "-") + "_" + fq2_list[1] + "_" + fq2_list[2] + "_" + fq2_list[3] + "_" + fq2_list[4] + ".fastq"
    #
    #     c1 = SeqIO.write(make_fastq_reads(fq1, row.forward_primer_sequence, row.reverse_primer_sequence), new_fq1, "fastq")
    #     print "Writing fastq file: " + new_fq1 + " (" + str(c1) + ")"
    #     c2 = SeqIO.write(make_fastq_reads(fq2, row.forward_primer_sequence, row.reverse_primer_sequence), new_fq2, "fastq")
    #     print "Writing fastq file: " + new_fq2 + " (" + str(c2) + ")"

if __name__ == "__main__":
    main(sys.argv)
