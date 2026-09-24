#!/usr/bin/env bash
# =============================================================================
# STEP 8 - Fig 5G / 5H : DUBR peaks vs. histone-mark ChIP-seq (deepTools).
#
# Fig 5G  "Heatmaps of ChIRP-seq and ChIP-seq signals at DUBR peaks grouped
#          into two clusters (C1 and C2) according to their signals by the
#          k-means algorithm."           -> C1 = 140 peaks, C2 = 1891 peaks
# Fig 5H  "Profile plots of ChIP-seq signal at DUBR peaks grouped in C1 and
#          C2 clusters" in WT (grey), dEX1 (red), dEX2 (blue).
#
# HOW IT WORKS
#   computeMatrix   : for every peak, take the signal in a +/-2.5 kb window from
#                     each bigWig (1 DUBR ChIRP + 4 histone marks, WT).
#   plotHeatmap     : --kmeans 2 partitions the peaks into 2 clusters. We tell
#                     it to cluster on the *histone* columns only
#                     (--clusterUsingSamples 2 3 4 5); otherwise the (much
#                     stronger) ChIRP signal drives the clustering and you get
#                     "strong vs weak DUBR peaks" instead of the paper's
#                     "peaks at regulatory elements vs the rest".
#   --outFileSortedRegions : writes the peaks with their cluster label, which
#                     we split into C1.bed / C2.bed and re-use for Fig 5H.
#
# WHY THE FIRST-PASS HEATMAP DID NOT MATCH THE PAPER
#   It only had EVEN / ODD / Input. Fig 5G is an *integration* figure: the
#   interesting axis is the histone marks (from GSE268500, step 07), not the
#   ChIRP signal itself.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 08_heatmaps_fig5GH.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_fig5GH
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/dubr/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
peak_dir="$out/peaks"
bw_dir="$out/bigwig"
avg="$out/tracks/avg"
fig_dir="$out/figures"
threads=16

peaks="$peak_dir/DUBR_peaks_final.bed"
mkdir -p "$fig_dir" "$peak_dir/clusters"

# ---- Fig 5G: k-means (k=2) heatmap ------------------------------------------
echo "=== Fig 5G: computeMatrix at $(wc -l < "$peaks") DUBR peaks (+/-2.5 kb) ==="
computeMatrix reference-point --referencePoint center \
    -S "$bw_dir/DUBR_CPM.bw" \
       "$avg/H3K4me3_WT.bw" "$avg/H3K4me1_WT.bw" "$avg/H3K27ac_WT.bw" "$avg/H3K27me3_WT.bw" \
    -R "$peaks" \
    -b 2500 -a 2500 --binSize 50 \
    --samplesLabel DUBR H3K4me3 H3K4me1 H3K27ac H3K27me3 \
    --missingDataAsZero \
    -o "$fig_dir/fig5G_matrix.gz" -p "$threads"

plotHeatmap -m "$fig_dir/fig5G_matrix.gz" \
    --kmeans 2 --clusterUsingSamples 2 3 4 5 \
    --sortUsingSamples 4 --sortUsing mean \
    --colorMap Reds \
    --zMin 0 --zMax 4 0.8 0.5 1.0 0.5 \
    --yMin 0 \
    --refPointLabel "peak" --xAxisLabel "distance (bp)" \
    --heatmapHeight 12 --heatmapWidth 3 \
    --outFileSortedRegions "$peak_dir/clusters/DUBR_peaks_kmeans2.bed" \
    -out "$fig_dir/fig5G_heatmap_kmeans2.png" --dpi 200

# deepTools writes the cluster label in the last column ("cluster_1" ...).
# Clusters are numbered by decreasing mean signal, so cluster_1 = the
# histone-marked cluster = the paper's C1.
for k in 1 2; do
    awk -v k="cluster_$k" 'BEGIN{OFS="\t"} $NF==k {print $1,$2,$3,$4}' \
        "$peak_dir/clusters/DUBR_peaks_kmeans2.bed" > "$peak_dir/clusters/C${k}.bed"
    echo "  C$k: $(wc -l < "$peak_dir/clusters/C${k}.bed") peaks   (paper: C1=140, C2=1891)"
done

# ---- Fig 5H: per-mark profiles at C1 / C2 in WT, dEX1, dEX2 ----------------
echo "=== Fig 5H: profile plots ==="
for mark in H3K4me3 H3K4me1 H3K27ac H3K27me3; do
    computeMatrix reference-point --referencePoint center \
        -S "$avg/${mark}_WT.bw" "$avg/${mark}_EX1.bw" "$avg/${mark}_EX2.bw" \
        -R "$peak_dir/clusters/C1.bed" "$peak_dir/clusters/C2.bed" \
        -b 2500 -a 2500 --binSize 50 \
        --samplesLabel WT dEX1 dEX2 \
        --missingDataAsZero \
        -o "$fig_dir/fig5H_${mark}_matrix.gz" -p "$threads"

    # --perGroup: one panel per region set (C1, C2), one line per sample
    plotProfile -m "$fig_dir/fig5H_${mark}_matrix.gz" \
        --perGroup --colors grey red blue \
        --plotTitle "$mark" --yAxisLabel "ChIP-seq (CPM)" \
        --refPointLabel "peak" --regionsLabel C1 C2 \
        --plotHeight 6 --plotWidth 6 \
        -out "$fig_dir/fig5H_${mark}_profile.png" --dpi 200 \
        --outFileNameData "$fig_dir/fig5H_${mark}_profile.tab"
done

# ---- Fig 5G, sanity version: does the ChIRP itself look the same in C1/C2? --
plotProfile -m "$fig_dir/fig5G_matrix.gz" --perGroup \
    -out "$fig_dir/fig5G_profiles_by_sample.png" --dpi 150 \
    --regionsLabel "all peaks" 2>/dev/null || true

echo "Fig 5G -> $fig_dir/fig5G_heatmap_kmeans2.png"
echo "Fig 5H -> $fig_dir/fig5H_<mark>_profile.png"
echo "clusters -> $peak_dir/clusters/C1.bed C2.bed"
