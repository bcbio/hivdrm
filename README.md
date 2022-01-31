# hivdrm
Detect HIV Drug Resitant Mutations using amplicon sequencing data

A production-ready re-implementation of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7699007/,
https://github.com/Wei-Shao/HIV-DRLink.

## Install

- Install conda (if not already present): https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
- Install mamba: `conda install mamba -n base -c conda-forge`
- Clone hivdrm: `git clone https://github.com/bcbio/hivdrm.git`
- `cd hivdrm`
- Create conda environment to run hivdrm: `mamba env create -n hivdrm_production --file environment.yml`
- add hivdrm to PATH: `export PATH=/path/to/hivdrm:$PATH` in .bashrc or .bash_profile

## Run

```bash
conda activate hivdrm_production
cd /path/project
hivdrm.py \
--barcodes barcodes.csv \
--reference reference.edited.fasta \
--threads 10 \
r1.fq.gz r2.fq.gz
conda deactivate
```
Some clusters/batch systems require `source activate/deactivate` instead of `conda activate/deactivate`.

Example of barcodes.csv:
```
Sample_ID,Primers,F-Linkers,R-Linkers
S01,F1/R1,CGCCTG,GCCATG
S02,F1/R2,CGCCTG,TACAAG
S03,F1/R3,CGCCTG,ATTCCG
S04,F1/R4,CGCCTG,TCGGGA
S05,F1/R5,CGCCTG,GAATGA
S06,F1/R6,CGCCTG,GCCTAA
S07,F2/R1,CGTGAT,GCCATG
S08,F2/R2,CGTGAT,TACAAG
S09,F2/R3,CGTGAT,ATTCCG
S10,F2/R4,CGTGAT,TCGGGA
S11,F2/R5,CGTGAT,GAATGA
S12,F2/R6,CGTGAT,GCCTAA
S13,F3/R1,CTGATC,GCCATG
```

Example of reference.fasta:
```
>EF602219.1 HIV-1 isolate 1779 from South Africa pol protein (pol) gene, partial cds edited
CGCCTGAATCCATATAACACTCCAATATTTGCCATAAAAAAGAAGGACAGTACTAAGTGGAGAAAATTAGTAGATTTCAGGGAACTTAATAAAAGAACTCAAGACTTTTGGGAAGTTCAATTAGGAATACCACATCCAGCAGGATTAAAAAAGAAAAAATCAGTGACAGTACTGGATGTGGGGGATGCATATTTTTCAGTTCCTTTAGATGAAGGCTTCAGAAAATATACTGCATTCACCATACCTAGTATAAACAATGAAACACCAGGGATTAGATATCAATATAATGTGCTCCCAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCCTTTAGAGCAAGAAATCCAGAAATAGTCATCTATCAATATATGGATGACTTGTATGTGGGATCTGACTTAGAAATAGGGCAACATAGAGCAAAAATAGAGGAATTAAGAGCACATTTATTAGGGTGGGGATTTACCACWCCAGACAAGAAACATCAGAAGGAACCCCCATTTCTTTGGATGGGGTACGAACTCCATCCTGACAAATGGACAGTNNNNNNNNNNCTAGCAGGATGACTTCGATACCCATGGC
```

Barcode processing steps for libraries with multiple samples and blastn step benefit a lot from multithreading. Consider running on a server/cluster with 10 threads/20G RAM.

## Output

- DRM.xlsx - 1st sheet - DRM stats, then one sheet per sample
- freq.xlsx - barcode stats, per sample all and top allele frequencies
- `_hivdrm_tmp` - intermediate files from all the steps

## Uninstall

- `conda remove --name hivdrm_production --all`
