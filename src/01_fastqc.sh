#!/usr/bin/env bash
# =============================================================================
# STEP 1 - Quality check the raw reads with FastQC.
#
# FastQC looks at each FASTQ file and reports read quality, adapter content,
# etc. It does not need the genome. We give it both R1 and R2 of every sample.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 01_fastqc.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_fastqc
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
fastq_dir="$out/fastq"
qc_dir="$out/fastqc_raw"
threads=8

mkdir -p "$qc_dir"

# Run FastQC on all 12 files (R1 + R2 for each of the 6 samples).
fastqc -t "$threads" -o "$qc_dir" \
    "$fastq_dir"/Input_rep1_R1.fastq.gz "$fastq_dir"/Input_rep1_R2.fastq.gz \
    "$fastq_dir"/Input_rep2_R1.fastq.gz "$fastq_dir"/Input_rep2_R2.fastq.gz \
    "$fastq_dir"/EVEN_rep1_R1.fastq.gz  "$fastq_dir"/EVEN_rep1_R2.fastq.gz \
    "$fastq_dir"/EVEN_rep2_R1.fastq.gz  "$fastq_dir"/EVEN_rep2_R2.fastq.gz \
    "$fastq_dir"/ODD_rep1_R1.fastq.gz   "$fastq_dir"/ODD_rep1_R2.fastq.gz \
    "$fastq_dir"/ODD_rep2_R1.fastq.gz   "$fastq_dir"/ODD_rep2_R2.fastq.gz

# MultiQC combines all the FastQC reports into one easy-to-read summary.
# multiqc "$qc_dir" -o "$out/multiqc_raw"

echo "FastQC reports are in: $qc_dir"
