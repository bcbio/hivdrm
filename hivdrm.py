#!/usr/bin/env python3

# usage nivdrm.py r1.fq.gz r2.fq.gz

import os
import re
import csv
import sys
import gzip
import shutil
import pandas as pd
import argparse
import textwrap
import subprocess
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Blast import NCBIXML

hivdrm_work_dir = "_hivdrm_tmp"
s5_prefix = "s5_demultiplex"
s7_prefix = "s7_blast_result"
s8_prefix = "s8_consensus"
s9_prefix = "s9_sierrapy"
s10_prefix = "s10_drlink"
samples = []

def touch(file_name):
    cf = Path(file_name)
    cf.touch(exist_ok = True)

def s1_rc_and_concatenate(r1, r2):
    """Reverse complement r2 and merge with r1 into a long amplicon"""
    print("Step1 - RC and concatenate")

    control_file = os.path.join(hivdrm_work_dir, "step1.done")
    result_file = os.path.realpath(os.path.join(hivdrm_work_dir, "step1.fq.gz"))
    if os.path.exists(control_file) and os.path.exists(result_file):
        return result_file

    i = 0
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
    touch(control_file)
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
    print("Step2 - quality filter")
    i = 0
    control_file = os.path.join(hivdrm_work_dir, "step2.done")
    result_file = os.path.join(hivdrm_work_dir, "step2.fq.gz")
    if os.path.exists(control_file) and os.path.exists(result_file):
        return result_file

    with gzip.open(file_fq_gz, "rt") as fq, gzip.open(result_file, "wt") as out:
        for record in SeqIO.parse(fq, "fastq"):
            if is_good_read(record, min_q = min_q, min_pct = min_pct):
                out.write(record.format("fastq"))
    touch(control_file)
    return result_file

def s3_fastq_to_fasta(file_fq_gz):
    """remove quality values"""
    print("Step3 - fastq to fasta ")
    control_file = os.path.join(hivdrm_work_dir, "step3.done")
    result_file = os.path.join(hivdrm_work_dir, "step3.fa")

    if os.path.exists(control_file) and os.path.exists(result_file):
        return result_file

    with gzip.open(file_fq_gz, "rt") as fq , open(result_file, "wt") as f_out:
        for record in SeqIO.parse(fq, "fastq"):
            SeqIO.write(record, f_out, "fasta-2line")
    touch(control_file)
    return result_file

def s4_trim_left_right(file_fa, bp_left = 4, bp_right = 4):
    """trim bp from left and right to align barcodes exactly to the edges"""
    print("Step4 - Trim left right ")
    control_file = os.path.join(hivdrm_work_dir, "step4.done")
    result_file = os.path.join(hivdrm_work_dir, "step4.fa")

    if os.path.exists(control_file) and os.path.exists(result_file):
        return result_file

    with open(file_fa, "rt") as fa_in, open(result_file, "wt") as fa_out:
        for record in SeqIO.parse(fa_in, "fasta-2line"):
            if len(record.seq) > bp_left + bp_right:
                record.seq = record.seq[bp_left:][:-bp_right]
                SeqIO.write(record, fa_out, "fasta-2line")
    touch(control_file)
    return result_file

def s5_demultiplex_samples(file_fa, barcodes_csv):
    """ demultiplex samples based on dual barcodes """
    print("Step5 - demultiplex samples")
    samples_barcode = {}
    result_file = "step5.stats.csv"
    control_file = os.path.join(hivdrm_work_dir, "step5.done")

    # still populate samples even the demultiplex work is done
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
            samples_barcode[combined_barcode] = sample_name
            samples.append(sample_name)

    os.makedirs(os.path.join(hivdrm_work_dir, s5_prefix), exist_ok = True)

    if os.path.exists(control_file) and os.path.exists(result_file):
        return result_file

    filedata = {filename: open(os.path.join(hivdrm_work_dir, s5_prefix, f"{filename}.fa"), "wt") for filename in samples_barcode.values()}
    unmatched = open(os.path.join(hivdrm_work_dir, s5_prefix, "unmatched.fa"), "wt")

    with open(file_fa, "rt") as fa_in:
        for record in SeqIO.parse(fa_in, "fasta-2line"):
            barcode = str(record.seq[0:barcode_f_len] + record.seq[-barcode_r_len:])
            if barcode in samples_barcode:
                SeqIO.write(record, filedata[samples_barcode[barcode]], "fasta-2line")
            else:
                SeqIO.write(record, unmatched, "fasta-2line")

    for f in filedata.values():
        f.close()
    unmatched.close()
    touch(control_file)
    touch(result_file)
    return result_file

def s6_create_blast_db(reference_fasta):
    print("Step6 - create blast db ")
    if shutil.which("makeblastdb") is None:
        print("makeblastdb is not found in the PATH. Load blast module or install NCBI blast")
        exit(1)
    blastdb_dir = os.path.join(hivdrm_work_dir, "blastdb")
    os.makedirs(blastdb_dir, exist_ok = True)
    bname = os.path.basename(reference_fasta)
    blastdb_path = os.path.join(blastdb_dir, bname)

    control_file = os.path.join(hivdrm_work_dir, "step6.done")
    if os.path.exists(control_file) and os.path.exists(blastdb_path):
        return blastdb_path

    shutil.copyfile(reference_fasta, blastdb_path)
    #prev_dir = os.getcwd()
    #os.chdir(blastdb_dir)
    cmd = f"makeblastdb -in {blastdb_path} -dbtype nucl"
    subprocess.check_call(cmd, shell = True)
    #os.chdir(prev_dir)
    touch(control_file)
    return blastdb_path

def s7_blastn_xml(qry, base):
    """ run blastn with xml output qry and base are absolute paths """
    print("Step7 ... :" + qry)
    os.makedirs(os.path.join(hivdrm_work_dir, s7_prefix), exist_ok = True)
    qry_file = os.path.basename(qry)
    base_file = os.path.basename(base)
    sample_name = os.path.splitext(qry_file)[0]
    result_file = f"{sample_name}.xml"
    result_path = os.path.realpath(os.path.join(hivdrm_work_dir, s7_prefix, result_file))
    if os.path.exists(result_path):
        return result_file
    cmd = (f"blastn -num_threads 1 "
           f"-query {qry} "
           f"-db {base} "
           f"-out {result_path} " \
           f"-dust no " \
           f"-num_alignments 1 " \
           f"-outfmt 5")
    subprocess.check_call(cmd, shell = True)
    return result_file

def s7_blast_all(barcodes_csv, reference_path):
    with open(barcodes_csv, "rt") as csvfile:
        csvreader = csv.reader(csvfile)
        # skip header
        next(csvfile, None)
        for row in csvreader:
            sample_fa = row[0] + ".fa"
            qry_path = os.path.realpath(os.path.join(hivdrm_work_dir, s5_prefix, sample_fa))
            s7_blastn_xml(qry_path, reference_path)

def hamming_dist(s1, s2):
    return sum(1 for (a, b) in zip(s1, s2) if a != b)

def get_umi(s_umi, dict_umi):
    """find 1-bp distant umi in a dictionary """
    s_result = s_umi
    if s_umi not in dict_umi:
        for s in dict_umi:
            dist = hamming_dist(s, s_umi)
            #sm = SequenceMatcher(None, s, s_umi)
            #if sm.ratio() >= 0.9:
            if dist == 1:
                s_result = s
                break
    return s_result

def get_consensus(family):
    s_consensus = ""
    for i in range(0, len(family[0])):
        s_column = ""
        for s in family:
            s_column += s[i]
        res = dict(Counter(s_column).most_common())
        top_allele = list(res.keys())[0]
        top_freq = res[top_allele]/len(family)
        #reference has an insert - but deleting makes fragments of different lengths
        #if top_allele == "-" and top_freq >= 0.5:
        #    continue
        n_alleles = len(list(res.keys()))
        if top_freq == 1.0 and n_alleles == 1:
            s_consensus += top_allele
        elif n_alleles > 1:
            second_allele = list(res.keys())[1]
            s_consensus += second_allele
        else:
            print(f"Error: {i} {s_column} {top_freq} {res} {family}")
            s_consensus += "N"
            # else:
            # s_consensus += "N"
    return s_consensus

def s8_make_consensus(sample, min_family_size = 5, max_family_size = 20):
    os.makedirs(os.path.join(hivdrm_work_dir, s8_prefix), exist_ok = True)
    sample_xml = sample + ".xml"
    blast_xml = os.path.realpath(os.path.join(hivdrm_work_dir, s7_prefix, sample_xml))
    result_handle = open(blast_xml)
    blast_records = NCBIXML.parse(result_handle)

    families = {}

    i=0
    for blast_record in blast_records:
        # taking one record one hsp
        try:
            hsp = blast_record.alignments[0].hsps[0]
            i += 1
        except:
            continue

        # hsp.align_length
        # ready query with gaps
        # hsp.query
        # reference with NNN
        # hsp.sbjct
        # UMI position
        m = re.search("N{5,}", hsp.sbjct)
        if not m:
            continue

        start = m.start()
        end = m.end() # end + 1 in fact
        # UMI is aligned with NNNNN
        raw_umi = hsp.query[start:end]
        umi = get_umi(raw_umi, families)
        #umi = raw_umi
            #print(f"HSP: {hsp.sbjct}")
                    #print(f"UMI: {umi}")
        seq = hsp.query[0:start]
        ref = hsp.sbjct[0:start]
        # we remove insertions in the query
        while ref.find("-") > -1:
            pos = ref.find("-")
            ref = ref[0:pos] + ref[pos+1:]
            seq = seq[0:pos] + seq[pos+1:]

        if umi in families:
            families[umi].append(seq)
        else:
            families[umi] = [ seq ]

        if i%1000 == 1:
            print(f"Total families: {len(families)}")
            print(f"Blast records processed: {i}")
        #if i > 10000:
        #    break
    results = {}
    consensus_fasta = sample + ".consensus.fasta"
    result_file = os.path.join(hivdrm_work_dir, s8_prefix, consensus_fasta)
    if os.path.exists(result_file):
        os.remove(result_file)

    for umi in families:
        filename = umi + ".fasta"
        family = families[umi]
        if len(family) >= min_family_size and len(family) <= max_family_size:
            lens = {len(s) for s in family}
            if len(lens) == 1:
                # lens are the same
                s_consensus = get_consensus(family)
                with open(result_file, "a") as f:
                    f.write(f">UMI_{umi}_FAMSIZE_{len(family)}_LEN_{len(s_consensus)}\n")
                    f.write(s_consensus + "\n")
        #fname = f"LEN_{len(s_consensus)}_{filename}"
        #with open(fname, "a") as f:
        #    for i, seq in enumerate(family):
        #        f.write(f">{umi}{i}\n")
        #        f.write(seq + "\n")
            # save a umi family for debug
            #else:
            #    with open(filename, "a") as f:
            #    for i, seq in enumerate(family):
            #        f.write(f">{umi}{i}\n")
            #        f.write(seq + "\n")
    result_handle.close()

def s8_make_consensus_all():
    control_file = os.path.join(hivdrm_work_dir, "step8.done")
    if os.path.exists(control_file):
        return
    for s in samples:
        s8_make_consensus(s)
    touch(control_file)

def s9_write_simple_mutations():
    os.makedirs(os.path.join(hivdrm_work_dir, s9_prefix), exist_ok = True)
    simple_mutations = os.path.join(hivdrm_work_dir, s9_prefix, "simple_mutations.gql")
    if os.path.exists(simple_mutations):
        return
    with open(simple_mutations, "wt") as fout: 
        fout.writelines(textwrap.dedent("""\
                           inputSequence {
                           header,
                           },
                           mutations {
                           gene {name}
                           primaryType
                           text
                           }"""))

def s9_sierrapy(sample):
    os.makedirs(os.path.join(hivdrm_work_dir, s9_prefix), exist_ok = True)
    consensus_fasta = sample + ".consensus.fasta"
    consensus_path = os.path.join(hivdrm_work_dir, s8_prefix, consensus_fasta)
    result_json = sample + ".json"
    result_path = os.path.join(hivdrm_work_dir, s9_prefix, result_json)
    if os.path.exists(result_path):
        return
    simple_mutations = os.path.join(hivdrm_work_dir, s9_prefix, "simple_mutations.gql")
    cmd = (f"sierrapy fasta -q {simple_mutations} "
           f"-o {result_path} {consensus_path}")
    subprocess.check_call(cmd, shell = True)

def s9_sierrapy_all():
    control_file = os.path.join(hivdrm_work_dir, "step9.done")
    if os.path.exists(control_file):
        return
    for sample in samples:
        s9_sierrapy(sample)
    touch(control_file)

def s10_drlink_perl(sample):
    print(f"step10: {sample}")
    s10_dir = os.path.join(hivdrm_work_dir, s10_prefix)
    os.makedirs(s10_dir, exist_ok = True)
    drlink_perl = os.path.join(os.path.dirname(__file__), "scripts", "HIV-DRLink_github.pl")
    result_json = sample + ".json"
    json_path = os.path.join(hivdrm_work_dir, s9_prefix, result_json)
    cmd = (f"perl {drlink_perl} {json_path}")
    drlink_stats = os.path.join(hivdrm_work_dir, "drlink_stats.tsv")
    with open(drlink_stats, "a") as fout:
        fout.write(sample + "\t")
        fout.seek(2)
        subprocess.check_call(cmd, shell = True, stdout = fout)
    result_file = sample + ".json_DRM.tsv"
    result_path = os.path.join(hivdrm_work_dir, s9_prefix, result_file)
    if os.path.exists(result_path):
        result_drm = sample + ".drm.tsv"
        shutil.move(result_path, os.path.join(s10_dir, result_drm))

# requires pip install openpyxl
def s10_drlink_all():
    control_file = os.path.join(hivdrm_work_dir, "step10.done")
    if os.path.exists(control_file):
        return
    for sample in samples:
        s10_drlink_perl(sample)
    # combine into one excel
    df = pd.DataFrame()
    s10_dir = os.path.join(hivdrm_work_dir, s10_prefix)
    tsvs = os.listdir(s10_dir)
    with pd.ExcelWriter("DRM.xlsx", engine = "openpyxl") as writer:
        tsv_path = os.path.join(hivdrm_work_dir, "drlink_stats.tsv")
        df = pd.read_csv(tsv_path, sep = '\t', names = ["sample", "summary"])
        df.to_excel(writer, sheet_name = "drlink_stats", index = False)
        for tsv_file in sorted(tsvs):
            if tsv_file.endswith(".tsv"):
                tsv_path = os.path.join(s10_dir, tsv_file)
                df = pd.read_csv(tsv_path, sep = '\t', skiprows = 1)
                df.to_excel(writer, sheet_name = tsv_file, index = False)
    touch(control_file)

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
    s8_make_consensus_all()
    s9_write_simple_mutations()
    s9_sierrapy_all()
    s10_drlink_all()

