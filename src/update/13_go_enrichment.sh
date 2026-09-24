#!/usr/bin/env bash
# =============================================================================
# STEP 13 - Slurm wrapper for 13_go_enrichment.R (GO enrichment: local
#           equivalent of the Metascape panels Fig 5F, 6F, S8A).
#
# Needs steps 09 and 11 (gene lists) and the chirp-fig environment with
# clusterProfiler (README section 2b).
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 13_go_enrichment.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_go
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/dubr/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chirp-fig

# Under sbatch, $0 is a spooled copy of this script -> use the submit folder.
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
Rscript 13_go_enrichment.R
