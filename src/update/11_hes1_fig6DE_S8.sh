#!/usr/bin/env bash
# =============================================================================
# STEP 11 - Fig 6D, 6E and Supp. Fig S8D, S8J (deepTools + bedtools).
#
#  Fig 6D  Heatmap + profiles of HES1, GATA-1 (ENCODE K562) and DUBR ChIRP at
#          the H3K27ac sites GAINED / REDUCED upon DUBR disruption (Table S4),
#          ranked by HES1 signal. Message of the panel: HES1 and GATA-1 sit on
#          these elements, DUBR itself does NOT (indirect effect).
#  Fig 6E  Donut: fraction of the 1216 shared DEGs that have a HES1 peak within
#          5 kb of their TSS (paper: 300 / 1216 = 25%).
#  Fig S8J Same as 6D with NCOR1 / RUNX1 / TAF1 (other direct DUBR targets that
#          are TFs).
#  Fig S8D ChIRP-seq profile at DUBR peaks linked to UP- vs DOWN-regulated
#          direct targets (BEDs written by step 09).
#
# Requires: step 07 (tracks), step 09 (Table S4 BEDs, DEG list, S8D BEDs).
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 11_hes1_fig6DE_S8.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_fig6
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
D="$out"
peak_dir="$D/peaks"; bw_dir="$D/bigwig"; fig_dir="$D/figures"
enc="$D/tracks/encode"; annot="$out/ref/annot"
threads=16

gained="$peak_dir/TableS4_H3K27ac_gained.bed"
reduced="$peak_dir/TableS4_H3K27ac_reduced.bed"
[ -s "$gained" ] || { echo "run 09_great_fig5CDE.R first (writes Table S4 BEDs)"; exit 1; }
echo "differential H3K27ac sites: gained $(wc -l < "$gained") | reduced $(wc -l < "$reduced")  (paper: 4000 | 1589)"

# ---- Fig 6D --------------------------------------------------------------------
echo "=== Fig 6D: HES1 / GATA-1 / DUBR at differential H3K27ac sites ==="
computeMatrix reference-point --referencePoint center \
    -S "$enc/HES1_K562_FC.bw" "$enc/GATA1_K562_FC.bw" "$bw_dir/DUBR_CPM.bw" \
    -R "$gained" "$reduced" \
    -b 4000 -a 4000 --binSize 100 \
    --samplesLabel HES1 GATA-1 DUBR --missingDataAsZero \
    -o "$fig_dir/fig6D_matrix.gz" -p "$threads"

# ENCODE tracks are fold-change-over-control (background ~1); DUBR is CPM.
# Leaving --zMax unset lets deepTools pick a per-sample 98th-percentile scale.
plotHeatmap -m "$fig_dir/fig6D_matrix.gz" \
    --sortUsingSamples 1 --sortUsing mean --sortRegions descend \
    --regionsLabel "Gained H3K27ac" "Reduced H3K27ac" \
    --colorMap Reds --zMin 0 --yMin 0 \
    --refPointLabel "peak" --xAxisLabel "distance (bp)" \
    --heatmapHeight 14 --heatmapWidth 3.5 --plotTitle "Fig 6D (ranked by HES1)" \
    -out "$fig_dir/fig6D_heatmap_HES1_GATA1_DUBR.png" --dpi 200

# ---- Fig S8J -------------------------------------------------------------------
echo "=== Fig S8J: NCOR1 / RUNX1 / TAF1 at differential H3K27ac sites ==="
computeMatrix reference-point --referencePoint center \
    -S "$enc/NCOR1_K562_FC.bw" "$enc/RUNX1_K562_FC.bw" "$enc/TAF1_K562_FC.bw" \
    -R "$gained" "$reduced" \
    -b 4000 -a 4000 --binSize 100 \
    --samplesLabel NCOR1 RUNX1 TAF1 --missingDataAsZero \
    -o "$fig_dir/figS8J_matrix.gz" -p "$threads"
plotHeatmap -m "$fig_dir/figS8J_matrix.gz" \
    --sortUsingSamples 1 --sortUsing mean --sortRegions descend \
    --regionsLabel "Gained H3K27ac" "Reduced H3K27ac" \
    --colorMap Oranges --zMin 0 --yMin 0 \
    --refPointLabel "peak" --heatmapHeight 14 --heatmapWidth 3.5 --plotTitle "Fig S8J (ranked by NCOR1)" \
    -out "$fig_dir/figS8J_heatmap_NCOR1_RUNX1_TAF1.png" --dpi 200

# ---- Fig 6E: DEGs with a HES1 peak within 5 kb of the TSS ----------------------
echo "=== Fig 6E: HES1-bound DEGs (HES1 ENCODE optimal IDR peaks, +/-5 kb of TSS) ==="
degs="$fig_dir/shared_DEGs_dEX1_dEX2.tsv"
zcat "$enc/HES1_K562_optimalIDR_peaks.bed.gz" | cut -f1-3 | sort -k1,1 -k2,2n > "$peak_dir/HES1_K562_peaks.bed"
echo "  HES1 peaks: $(wc -l < "$peak_dir/HES1_K562_peaks.bed")   (paper Fig S8H: 10,634)"
# TSS (one per transcript, UCSC refGene) of the shared DEGs
awk 'NR>1{print $1}' "$degs" | sort -u > "$fig_dir/tmp_deg_ids.txt"
awk 'NR==FNR{d[$1]=1; next} ($4 in d)' "$fig_dir/tmp_deg_ids.txt" "$annot/tss.bed" \
    | sort -k1,1 -k2,2n > "$peak_dir/shared_DEGs_TSS.bed"
n_deg=$(cut -f4 "$peak_dir/shared_DEGs_TSS.bed" | sort -u | wc -l)
bedtools window -w 5000 -u -a "$peak_dir/shared_DEGs_TSS.bed" -b "$peak_dir/HES1_K562_peaks.bed" \
    | cut -f4 | sort -u > "$fig_dir/fig6E_HES1_bound_DEGs.txt"
n_bound=$(wc -l < "$fig_dir/fig6E_HES1_bound_DEGs.txt")
echo "  DEGs with TSS annotation: $n_deg ; HES1-bound (<=5 kb): $n_bound   (paper: 1216 -> 300, 25%)"
rm -f "$fig_dir/tmp_deg_ids.txt"

python - "$n_deg" "$n_bound" "$fig_dir/fig6E_donut_HES1_bound_DEGs.png" <<'EOF'
import sys, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
n, k, outpng = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
fig, ax = plt.subplots(figsize=(4.5, 4.5))
ax.pie([k, n - k], colors=["#e08a3c", "#3b6fb6"], startangle=90, counterclock=False,
       wedgeprops=dict(width=0.35, edgecolor="white"))
ax.text(0, 0, f"dDUBR DEGs\n({n:,})", ha="center", va="center", fontsize=11)
ax.set_title(f"Fig 6E  HES1-bound genes: {k:,} ({100*k/n:.0f}%)\npaper: 300 / 1,216 (25%)", fontsize=10)
plt.tight_layout(); plt.savefig(outpng, dpi=200)
EOF
echo "  -> gene list for Metascape (Fig 6F): $fig_dir/fig6E_HES1_bound_DEGs.txt"

# ---- Fig S8D: ChIRP signal at DUBR peaks of up- vs down-regulated targets -------
up="$peak_dir/DUBR_peaks_direct_upregulated.bed"; dn="$peak_dir/DUBR_peaks_direct_downregulated.bed"
if [ -s "$up" ] && [ -s "$dn" ]; then
    echo "=== Fig S8D: ChIRP profile, peaks of up ($(wc -l < "$up")) vs down ($(wc -l < "$dn")) direct targets ==="
    computeMatrix reference-point --referencePoint center \
        -S "$bw_dir/DUBR_CPM.bw" -R "$up" "$dn" \
        -b 2500 -a 2500 --binSize 50 --samplesLabel "DUBR ChIRP" --missingDataAsZero \
        -o "$fig_dir/figS8D_matrix.gz" -p "$threads"
    plotProfile -m "$fig_dir/figS8D_matrix.gz" --colors "#2e4a9e" "#d62728" \
        --regionsLabel "Upregulated" "Downregulated" --yAxisLabel "ChIRP-seq (CPM)" \
        --refPointLabel "peak" --plotTitle "Fig S8D" --plotHeight 6 --plotWidth 6 \
        -out "$fig_dir/figS8D_chirp_profile_up_vs_down.png" --dpi 200
fi

echo "Fig 6D -> $fig_dir/fig6D_heatmap_HES1_GATA1_DUBR.png"
echo "Fig 6E -> $fig_dir/fig6E_donut_HES1_bound_DEGs.png"
echo "Fig S8J -> $fig_dir/figS8J_heatmap_NCOR1_RUNX1_TAF1.png ; S8D -> figS8D_chirp_profile_up_vs_down.png"
