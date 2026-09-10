#!/usr/bin/env bash
# =============================================================================
# STEP 6 - Make signal tracks and plots with deepTools (the final figures).
#
# This step does three things:
#   a) Turn each merged BAM into a "bigWig" signal track (a genome-browser file
#      showing read density), normalized so samples are comparable.
#   b) Check that EVEN and ODD signals agree (a correlation scatterplot).
#   c) Draw a heatmap and a profile plot of the signal around the DUBR peaks.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 06_signal_deeptools.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_deeptools
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/dubr/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
bam_dir="$out/bam"
peak_dir="$out/peaks"
bw_dir="$out/bigwig"       # output: signal tracks
fig_dir="$out/figures"     # output: plots
threads=16

mkdir -p "$bw_dir" "$fig_dir"

# ---- a) Make a normalized signal track (bigWig) for each condition ----------
# --normalizeUsing CPM makes samples comparable despite different read counts.
for cond in Input EVEN ODD; do
    bamCoverage \
        -b "$bam_dir/${cond}_merged.bam" \
        -o "$bw_dir/${cond}_CPM.bw" \
        --normalizeUsing CPM \
        --binSize 10 \
        -p "$threads"
done

# ---- b) Do EVEN and ODD agree? Measure signal at the peaks, then correlate --
multiBigwigSummary BED-file \
    --BED "$peak_dir/DUBR_confident_peaks.bed" \
    -b "$bw_dir/EVEN_CPM.bw" "$bw_dir/ODD_CPM.bw" \
    --labels EVEN ODD \
    -o "$fig_dir/EVEN_ODD_peakcounts.npz" \
    --outRawCounts "$fig_dir/EVEN_ODD_peakcounts.tab" \
    -p "$threads"

plotCorrelation \
    -in "$fig_dir/EVEN_ODD_peakcounts.npz" \
    --corMethod pearson --whatToPlot scatterplot \
    --labels EVEN ODD \
    -o "$fig_dir/EVEN_ODD_correlation.png"

# ---- c) Heatmap + profile of signal around the confident peaks --------------
# computeMatrix gathers the signal in a window 3000 bp on each side of each peak
# center; plotHeatmap and plotProfile then draw it.
computeMatrix reference-point \
    --referencePoint center \
    -S "$bw_dir/EVEN_CPM.bw" "$bw_dir/ODD_CPM.bw" "$bw_dir/Input_CPM.bw" \
    -R "$peak_dir/DUBR_confident_peaks.bed" \
    -a 3000 -b 3000 \
    --skipZeros \
    -o "$fig_dir/DUBR_peaks_matrix.gz" \
    -p "$threads"

plotHeatmap \
    -m "$fig_dir/DUBR_peaks_matrix.gz" \
    --kmeans 2 \
    --colorMap Reds \
    --refPointLabel "peak center" \
    -out "$fig_dir/DUBR_peaks_heatmap.png"

plotProfile \
    -m "$fig_dir/DUBR_peaks_matrix.gz" \
    -out "$fig_dir/DUBR_peaks_profile.png"

echo "Signal tracks are in: $bw_dir"
echo "Plots are in:         $fig_dir"
