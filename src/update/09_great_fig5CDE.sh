#!/usr/bin/env bash
# =============================================================================
# STEP 9 - Slurm wrapper for the R script 09_great_fig5CDE.R
#          (Fig 5C, 5D, 5E, EVEN/ODD scatter, Fig S8C, Table S4 BEDs for step 11).
#
# Slurm can only submit shell scripts, so this file just turns on the R
# environment and runs the .R file that sits next to it.
#
# BEFORE RUNNING (once): copy the paper's supplementary tables to scratch:
#     mkdir -p /tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr/paper
#     scp doc/dubr_paper/Table_S1-S6.xlsx \
#         <you>@login.tscc.sdsc.edu:/tscc/lustre/ddn/scratch/<you>/chirp_seq_analysis/dubr/paper/
#
# NOTE: the first run downloads GREAT's hg19 TSS table (a few MB) from GitHub,
# so the node needs internet. If that fails, run on the login node:
#     conda activate chirp-fig && Rscript 09_great_fig5CDE.R
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 09_great_fig5CDE.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_great
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/dubr/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chirp-fig          # R + rGREAT + ComplexHeatmap ... (README section 2b)

# Under sbatch, $0 points to a spooled copy of this script, so use the folder
# the job was submitted from (SLURM_SUBMIT_DIR) to find the .R file.
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
Rscript 09_great_fig5CDE.R
