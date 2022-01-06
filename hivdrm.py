#!/usr/bin/env python3

# usage nivdrm.py r1.fq.gz r2.fq.gz

import sys
import gzip
from Bio.Seq import Seq
from Bio import SeqIO

def s1_rc_and_concatenate(r1, r2):
    """Reverse complement r2 and merge with r1 into a long amplicon"""
    i = 0
    result_file = "step1.fq.gz"
    with gzip.open(sys.argv[1], "rt") as fq1, \
         gzip.open(sys.argv[2], "rt") as fq2, \
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

if __name__ == "__main__":
    description = "Detect HIV drug resistant mutations"
    r1 = sys.argv[1]
    r2 = sys.argv[2]
    step1_out = s1_rc_and_concatenate(r1, r2)
    step2_out = s2_fastq_quality_filter(step1_out)
