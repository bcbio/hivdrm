#!/usr/bin/env python3

# usage nivdrm.py r1.fq.gz r2.fq.gz

import os
import csv
import sys
import gzip
import argparse

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO

def s1_rc_and_concatenate(r1, r2):
    """Reverse complement r2 and merge with r1 into a long amplicon"""
    i = 0
    result_file = "step1.fq.gz"
    with gzip.open(r1, "rt") as fq1, \
         gzip.open(r2, "rt") as fq2, \
         gzip.open(result_file, "wt") as out:
        for line1, line2 in zip(fq1, fq2):
            line1 = line1.strip()
            line2 = line2.strip()
            out_buf = ""
            if i % 4 == 0 or i % 4 == 2 :
                out_buf = line1
            elif i % 4 == 1:
                seq2 = Seq(line2)
                out_buf = line1 + str(seq2.reverse_complement())
            elif i % 4 == 3:
                out_buf = line1 + line2[::-1]
                i = -1
            out.write(out_buf + "\n")
            i += 1
    return result_file

def is_good_read(a_read, min_q, min_pct):
    """Check whether a read has >= min_pct bases of Q >= min_q"""
    l_qual = a_read.letter_annotations["phred_quality"]
    good_bases = list(filter(lambda q: q >= min_q, l_qual))
    good_pct = round(100 * len(good_bases) / len(l_qual))
    result = good_pct >= min_pct
    return result

def s2_fastq_quality_filter(file_fq_gz, min_q = 20, min_pct = 90):
    """implementation of fastx_toolkit/fastq_quality_filter
       file is of single reads - amplicons  """
    i = 0
    result_file = "step2.fq.gz"
    with gzip.open(file_fq_gz, "rt") as fq, gzip.open(result_file, "wt") as out:
        for record in SeqIO.parse(fq, "fastq"):
            if is_good_read(record, min_q = min_q, min_pct = min_pct):
                out.write(record.format("fastq"))
    return result_file

def s3_fastq_to_fasta(file_fq_gz):
    """remove quality values"""
    result_file = "step3.fa"
    with gzip.open(file_fq_gz, "rt") as fq , open(result_file, "wt") as f_out:
        for record in SeqIO.parse(fq, "fastq"):
            SeqIO.write(record, f_out, "fasta-2line")
    return result_file

def s4_trim_left_right(file_fa, bp_left = 4, bp_right = 4):
    """trim bp from left and right to align barcodes exactly to the edges"""
    result_file = "step4.fa"
    with open(file_fa, "rt") as fa_in, open(result_file, "wt") as fa_out:
        for record in SeqIO.parse(fa_in, "fasta-2line"):
            if len(record.seq) > bp_left + bp_right:
                record.seq = record.seq[bp_left:][:-bp_right]
                SeqIO.write(record, fa_out, "fasta-2line")
    return result_file

def s5_demultiplex_samples(file_fa, barcodes_csv):
    """ demultiplex samples based on dual barcodes """
    samples = {}
    result_file = "step5.stats.csv"
    with open(barcodes_csv) as csvfile:
        csvreader = csv.reader(csvfile)
        n_row = len(csvfile.readlines())
        csvfile.seek(0)
        header = list(next(csvreader, None))
        n_col = len(header)
        if n_col != 4 or n_row == 0:
            print("Expecting --barcodes barcodes.csv to be a 4 column csv file: sampleid,primer_id,f-linker,r-linker")
            exit(1)
        else:
            header = list(next(csvreader, None))
            barcode_f_len = len(header[2])
            barcode_r_len = len(header[3])
        csvfile.seek(0)
        # skip header
        next(csvfile, None)
        for row in csvreader:
            sample_name = row[0]
            f_barcode = row[2]
            r_barcode = row[3]
            s = Seq(r_barcode)
            combined_barcode = f_barcode + s.reverse_complement()
            samples[combined_barcode] = sample_name

    s5_prefix = "s5_demultiplex"
    os.makedirs(s5_prefix, exist_ok = True)

    filedata = {filename: open(os.path.join(s5_prefix, f"{filename}.fa"), "wt") for filename in samples.values()}
    unmatched = open(os.path.join(s5_prefix, "unmatched.fa"), "wt")

    with open(file_fa, "rt") as fa_in:
        for record in SeqIO.parse(fa_in, "fasta-2line"):
            barcode = str(record.seq[0:barcode_f_len] + record.seq[-barcode_r_len:])
            if barcode in samples:
                SeqIO.write(record, filedata[samples[barcode]], "fasta-2line")
            else:
                SeqIO.write(record, unmatched, "fasta-2line")

    for f in filedata.values():
        f.close()
    unmatched.close()
    
    return result_file

def get_args(description):
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument("--barcodes", required = True)
    parser.add_argument("fastq_files", nargs = 2)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    description = "Detect HIV drug resistant mutations"
    args = get_args(description)
    r1 = args.fastq_files[0]
    r2 = args.fastq_files[1]
    step1_out = s1_rc_and_concatenate(r1, r2)
    step2_out = s2_fastq_quality_filter(step1_out)
    step3_out = s3_fastq_to_fasta(step2_out)
    step4_out = s4_trim_left_right(step3_out)
    step5_out = s5_demultiplex_samples(step4_out, args.barcodes)
