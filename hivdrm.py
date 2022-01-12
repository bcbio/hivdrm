#!/usr/bin/env python3

# usage nivdrm.py r1.fq.gz r2.fq.gz

import os
import csv
import sys
import gzip
import shutil
import argparse
import subprocess

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO

hivdrm_work_dir = "_hivdrm_tmp"
s5_prefix = "s5_demultiplex"

def s1_rc_and_concatenate(r1, r2):
    """Reverse complement r2 and merge with r1 into a long amplicon"""
    i = 0
    result_file = os.path.join(hivdrm_work_dir, "step1.fq.gz")
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
    result_file = os.path.join(hivdrm_work_dir, "step2.fq.gz")
    with gzip.open(file_fq_gz, "rt") as fq, gzip.open(result_file, "wt") as out:
        for record in SeqIO.parse(fq, "fastq"):
            if is_good_read(record, min_q = min_q, min_pct = min_pct):
                out.write(record.format("fastq"))
    return result_file

def s3_fastq_to_fasta(file_fq_gz):
    """remove quality values"""
    result_file = os.path.join(hivdrm_work_dir, "step3.fa")
    with gzip.open(file_fq_gz, "rt") as fq , open(result_file, "wt") as f_out:
        for record in SeqIO.parse(fq, "fastq"):
            SeqIO.write(record, f_out, "fasta-2line")
    return result_file

def s4_trim_left_right(file_fa, bp_left = 4, bp_right = 4):
    """trim bp from left and right to align barcodes exactly to the edges"""
    result_file = os.path.join(hivdrm_work_dir, "step4.fa")
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
    with open(barcodes_csv, "rt") as csvfile:
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

    os.makedirs(os.path.join(hivdrm_work_dir, s5_prefix), exist_ok = True)

    filedata = {filename: open(os.path.join(hivdrm_work_dir, s5_prefix, f"{filename}.fa"), "wt") for filename in samples.values()}
    unmatched = open(os.path.join(hivdrm_work_dir, s5_prefix, "unmatched.fa"), "wt")

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

def s6_create_blast_db(reference_fasta):
    if shutil.which("makeblastdb") is None:
        print("makeblastdb is not found in the PATH. Load blast module or install NCBI blast")
        exit(1)
    blastdb_dir = os.path.join(hivdrm_work_dir, "blastdb")
    os.makedirs(blastdb_dir, exist_ok = True)
    bname = os.path.basename(reference_fasta)
    blastdb_path = os.path.join(blastdb_dir, bname)
    shutil.copyfile(reference_fasta, blastdb_path)
    #prev_dir = os.getcwd()
    #os.chdir(blastdb_dir)
    cmd = f"makeblastdb -in {blastdb_path} -dbtype nucl"
    subprocess.check_call(cmd, shell = True)
    #os.chdir(prev_dir)
    return blastdb_path

def s7_blastn_xml(qry, base):
    """ run blastn with xml output qry and base are absolute paths """
    s7_prefix = "s7_blast_result"
    os.makedirs(os.path.join(hivdrm_work_dir, s7_prefix), exist_ok = True)
    qry_file = os.path.basename(qry)
    base_file = os.path.basename(base)
    result_file = f"{qry_file}_vs_{base_file}.xml"
    result_path = os.path.realpath(os.path.join(hivdrm_work_dir, s7_prefix, result_file))
    if os.path.isfile(result_path):
        os.remove(result_path)
    cmd = (f"blastn -num_threads 1 "
           f"-query {qry} "
           f"-db {base} "
           f"-out {result_path} " \
           f"-dust no " \
           f"-num_alignments 1 " \
           f"-outfmt 5")
    subprocess.check_call(cmd, shell = True)

def s7_blast_all(barcodes_csv, reference_path):
    with open(barcodes_csv, "rt") as csvfile:
        csvreader = csv.reader(csvfile)
        # skip header
        next(csvfile, None)
        for row in csvreader:
            sample_fa = row[0] + ".fa"
            qry_path = os.path.realpath(os.path.join(hivdrm_work_dir, s5_prefix, sample_fa))
            s7_blastn_xml(qry_path, reference_path)

def get_args(description):
    parser = argparse.ArgumentParser(description = description, usage = "%(prog)s [options]")
    parser.add_argument("--barcodes", required = True, help = "barcodes.csv")
    parser.add_argument("--reference", required = True, help = "reference.fasta")
    parser.add_argument("fastq_files", nargs = 2, help = "f1.fq.gz f2.fq.gz")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    description = "Detect HIV drug resistant mutations"
    args = get_args(description)
    r1 = args.fastq_files[0]
    r2 = args.fastq_files[1]
    os.makedirs(hivdrm_work_dir, exist_ok = True)
    step1_out = s1_rc_and_concatenate(r1, r2)
    step2_out = s2_fastq_quality_filter(step1_out)
    step3_out = s3_fastq_to_fasta(step2_out)
    step4_out = s4_trim_left_right(step3_out)
    step5_out = s5_demultiplex_samples(step4_out, args.barcodes)
    s6_out_fasta_ref = s6_create_blast_db(args.reference)
    s7_blast_all(args.barcodes, s6_out_fasta_ref)

