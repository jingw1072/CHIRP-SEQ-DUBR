#!/usr/bin/env bash
# =============================================================================
# STEP 1b - Clean up the reads with fastp (trim adapters and poly-G tails).
#
# The raw FastQC showed leftover sequencing adapters and long runs of "G"
# (a common artifact of Illumina machines). fastp removes both so they don't
# confuse the alignment step. Because these reads are short (38 bp), we throw
# away anything shorter than 20 bp after trimming.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 01b_trim_fastp.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_fastp
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/dubr/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu
# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
fastq_dir="$out/fastq"           # input: raw reads from step 00
trim_dir="$out/fastq_trimmed"    # output: cleaned reads (used by step 02)
rep_dir="$out/fastp_reports"     # fastp's html/json reports and logs
threads=8

mkdir -p "$trim_dir" "$rep_dir"

# Loop over each sample and clean up its two read files.
for sample in Input_rep1 Input_rep2 EVEN_rep1 EVEN_rep2 ODD_rep1 ODD_rep2; do
    echo "=== trimming $sample ==="
    fastp \
        --in1  "$fastq_dir/${sample}_R1.fastq.gz" \
        --in2  "$fastq_dir/${sample}_R2.fastq.gz" \
        --out1 "$trim_dir/${sample}_R1.fastq.gz" \
        --out2 "$trim_dir/${sample}_R2.fastq.gz" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --length_required 20 \
        --thread "$threads" \
        --html "$rep_dir/${sample}_fastp.html" \
        --json "$rep_dir/${sample}_fastp.json" \
        2> "$rep_dir/${sample}_fastp.log"
    echo "Done: $sample"
done

echo "Cleaned reads are in: $trim_dir"
