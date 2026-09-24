#!/usr/bin/env bash
# =============================================================================
# STEP 5 - Find the high-confidence DUBR peaks (peaks found in BOTH EVEN & ODD).
#
# EVEN and ODD are two independent sets of probes for the same DUBR RNA. A real
# binding site should show up in BOTH. So we keep only the peaks that appear in
# the EVEN set AND the ODD set (overlapping by more than 100 bp). These are our
# final, trustworthy DUBR peaks.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 05_confident_peaks.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_confident
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/dubr/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
peak_dir="$out/peaks"
fc_min=5
min_overlap=100     # the two peaks must overlap by more than 100 bp to count

even_np="$peak_dir/DUBR_EVEN_FC${fc_min}_noBL.narrowPeak"
odd_np="$peak_dir/DUBR_ODD_FC${fc_min}_noBL.narrowPeak"

# Trim each peak file down to 6 columns (chromosome, start, end, name,
# fold-enrichment, q-value) and sort by position. The 6-column format makes the
# overlap width land in a predictable column (13) after the intersect below.
#-k1,1 means sort using column 1 only (usually chromosome). 
#-k2,2n means within each chromosome, sort by column 2 numerically (n = numeric)
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$7,$9}' "$even_np" \
    | sort -k1,1 -k2,2n > "$peak_dir/EVEN_6col.bed"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$7,$9}' "$odd_np" \
    | sort -k1,1 -k2,2n > "$peak_dir/ODD_6col.bed"

# Find EVEN peaks that overlap ODD peaks. "-wo" also reports the overlap width
# (in the last column, 13). We then keep only overlaps wider than 100 bp.
bedtools intersect -wo -a "$peak_dir/EVEN_6col.bed" -b "$peak_dir/ODD_6col.bed" \
    | awk -v m="$min_overlap" 'BEGIN{OFS="\t"} $13 > m' \
    > "$peak_dir/Intersect_EVENvsODD_FC5-FDR0.01.bed"

# Keep just the peak location (first 4 columns), remove duplicates -> final set.
cut -f1-4 "$peak_dir/Intersect_EVENvsODD_FC5-FDR0.01.bed" \
    | sort -k1,1 -k2,2n -u \
    > "$peak_dir/DUBR_confident_peaks.bed"

# ---- Numbers for the Fig 5C Venn: EVEN-only / common / ODD-only --------------
# "common" is counted on each side separately (one EVEN peak can overlap two
# ODD peaks and vice versa), so EVEN-only = EVEN total - EVEN peaks with a partner.
n_even=$(wc -l < "$even_np")
n_odd=$(wc -l < "$odd_np")
n_even_common=$(cut -f4 "$peak_dir/Intersect_EVENvsODD_FC5-FDR0.01.bed" | sort -u | wc -l)
n_odd_common=$(cut -f10 "$peak_dir/Intersect_EVENvsODD_FC5-FDR0.01.bed" | sort -u | wc -l)
printf "set\tcount\nEVEN_total\t%d\nODD_total\t%d\nEVEN_only\t%d\nODD_only\t%d\ncommon_EVEN_side\t%d\ncommon_ODD_side\t%d\n" \
    "$n_even" "$n_odd" $((n_even - n_even_common)) $((n_odd - n_odd_common)) \
    "$n_even_common" "$n_odd_common" > "$peak_dir/fig5C_venn_counts.tsv"

# Report the counts at each stage.
echo "EVEN peaks:            $n_even"
echo "ODD peaks:             $n_odd"
echo "Overlapping (>100 bp): $(wc -l < "$peak_dir/Intersect_EVENvsODD_FC5-FDR0.01.bed")"
echo "Confident DUBR peaks:  $(wc -l < "$peak_dir/DUBR_confident_peaks.bed")"
echo "Fig 5C: EVEN-only $((n_even - n_even_common)) | common $n_even_common | ODD-only $((n_odd - n_odd_common))"
echo "   (paper Fig. 5C: EVEN-only 783 | common 2031 | ODD-only 1253)"

echo "Final peaks are in: $peak_dir/DUBR_confident_peaks.bed"
