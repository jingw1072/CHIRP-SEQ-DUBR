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
# We use the UCSC-named BAMs from step 03b so the tracks say "chr1" like every
# public track we will overlay later.
for cond in Input EVEN ODD; do
    bamCoverage \
        -b "$bam_dir/${cond}_merged.chr.bam" \
        -o "$bw_dir/${cond}_CPM.bw" \
        --normalizeUsing CPM \
        --binSize 10 \
        -p "$threads"
done
# One combined "DUBR" track = average of EVEN and ODD (both target the same RNA).
bigwigAverage -b "$bw_dir/EVEN_CPM.bw" "$bw_dir/ODD_CPM.bw" \
    -o "$bw_dir/DUBR_CPM.bw" --binSize 10 -p "$threads"

# ---- b) Do EVEN and ODD agree on HOW MUCH signal? -> final peak set ---------
# The paper keeps only peaks whose EVEN/ODD signal ratio is between 0.5 and 2
# (a peak seen strongly by one probe set but weakly by the other is suspect).
#   (1) multiBigwigSummary : mean CPM of EVEN / ODD / Input over each peak
#   (2) bedtools intersect : put the peak names back (column 4 is dropped in (1))
#   (3) awk                : ratio = EVEN/ODD, keep 0.5 <= ratio <= 2
conf="$peak_dir/DUBR_confident_peaks.bed"

multiBigwigSummary BED-file \
    --BED "$conf" \
    -b "$bw_dir/EVEN_CPM.bw" "$bw_dir/ODD_CPM.bw" "$bw_dir/Input_CPM.bw" \
    --labels EVEN ODD Input \
    -o "$peak_dir/EVEN_ODD_peak_signal.npz" \
    --outRawCounts "$peak_dir/EVEN_ODD_peak_signal.tab" \
    -p "$threads"

# The table has a '#' header line and may contain "nan" (awk: "nan"+0 = 0).
grep -v '^#' "$peak_dir/EVEN_ODD_peak_signal.tab" | sort -k1,1 -k2,2n \
    | bedtools intersect -wa -wb -f 1.0 -r -a "$conf" -b - \
    | awk 'BEGIN{OFS="\t"; print "chr","start","end","name","EVEN_cpm","ODD_cpm","Input_cpm","EVEN_over_ODD","pass_0.5_2"}
           { even=$8+0; odd=$9+0; inp=$10+0
             ratio = (odd > 0) ? even/odd : "inf"
             pass  = (odd > 0 && ratio >= 0.5 && ratio <= 2) ? 1 : 0
             print $1,$2,$3,$4,even,odd,inp,ratio,pass }' \
    > "$peak_dir/DUBR_peaks_EVEN_ODD_signal.tsv"

awk 'BEGIN{OFS="\t"} NR>1 && $9==1 {print $1,$2,$3,$4}' "$peak_dir/DUBR_peaks_EVEN_ODD_signal.tsv" \
    > "$peak_dir/DUBR_peaks_final.bed"

echo "  confident peaks       : $(wc -l < "$conf")"
echo "  pass 0.5<=EVEN/ODD<=2 : $(wc -l < "$peak_dir/DUBR_peaks_final.bed")  -> DUBR_peaks_final.bed"

# Pearson correlation of EVEN vs ODD on the final peaks (the paper's QC).
multiBigwigSummary BED-file \
    --BED "$peak_dir/DUBR_peaks_final.bed" \
    -b "$bw_dir/EVEN_CPM.bw" "$bw_dir/ODD_CPM.bw" --labels EVEN ODD \
    -o "$peak_dir/EVEN_ODD_final_signal.npz" -p "$threads"
plotCorrelation -in "$peak_dir/EVEN_ODD_final_signal.npz" \
    --corMethod pearson --whatToPlot scatterplot --log1p \
    --outFileCorMatrix "$fig_dir/QC_EVEN_ODD_pearson.txt" \
    -o "$fig_dir/QC_EVEN_ODD_scatter.png"
echo "  Pearson r (EVEN vs ODD, final peaks): $(grep "^'EVEN'" "$fig_dir/QC_EVEN_ODD_pearson.txt" | cut -f3)"

# ---- c) QC heatmap of EVEN / ODD / Input around the final peaks -------------
# (Fig 5G in the paper also needs the histone ChIP-seq tracks - see the
#  scripts 07-11 in src/dubr/ for the full figure reproduction.)
computeMatrix reference-point \
    --referencePoint center \
    -S "$bw_dir/EVEN_CPM.bw" "$bw_dir/ODD_CPM.bw" "$bw_dir/Input_CPM.bw" \
    -R "$peak_dir/DUBR_peaks_final.bed" \
    -a 3000 -b 3000 --binSize 50 \
    --samplesLabel EVEN ODD Input \
    -o "$fig_dir/QC_chirp_matrix.gz" \
    -p "$threads"

plotHeatmap \
    -m "$fig_dir/QC_chirp_matrix.gz" \
    --sortUsingSamples 1 2 --sortUsing mean \
    --colorMap Reds --zMax 3 --yMax 4 \
    --refPointLabel "peak center" \
    -out "$fig_dir/QC_chirp_EVEN_ODD_Input_heatmap.png"

echo "Signal tracks are in: $bw_dir"
echo "Final peak set:       $peak_dir/DUBR_peaks_final.bed"
echo "Plots are in:         $fig_dir"
