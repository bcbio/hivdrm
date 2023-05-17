#!/bin/bash
# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=24:00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=hivdrm          # Job name
#SBATCH -c 20			    # cores
#SBATCH --mem=20G               # Memory needed per CPU or --mem-per-cpu
#SBATCH --output=hivdrm_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=hivdrm_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

date

. ~/work/bcbio_devel/.bcbio_devel_profile

source activate hivdrm_production

~/code/hivdrm/hivdrm.py \
--barcodes $1 \
--reference $2 \
--threads 20 \
$3 $4

source deactivate

date
