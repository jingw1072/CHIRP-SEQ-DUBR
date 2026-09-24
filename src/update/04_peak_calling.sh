#!/usr/bin/env bash
# =============================================================================
# STEP 4 - Call peaks with MACS2 (find where DUBR binds the genome).
#
# A "peak" is a spot where the EVEN or ODD sample has many more reads than the
# Input (background) sample - i.e. a likely DUBR binding site. We then keep only
# strong peaks (at least 5x more signal than background) and remove peaks that
# fall in known problem regions (the blacklist).
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 04_peak_calling.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_macs2
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
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
blacklist="$out/ref/hg19-blacklist.v2.bed"
fc_min=5                          # keep peaks with >= 5x signal over background

mkdir -p "$peak_dir"

# Do this for both probe pools: EVEN and ODD, each compared to Input.
for pool in EVEN ODD; do
    echo "=== MACS2 callpeak: $pool vs Input ==="
    # -t treatment (EVEN or ODD), -c control (Input)
    # -f BAMPE  the data is paired-end
    # -g hs     human genome
    # -q 0.01   statistical cutoff (false discovery rate < 1%)
    # Uses the UCSC-named BAMs from step 03b (chr1, chr2...) so that the
    # blacklist filter below actually matches chromosome names.
    macs2 callpeak \
        -t "$bam_dir/${pool}_merged.chr.bam" \
        -c "$bam_dir/Input_merged.chr.bam" \
        -f BAMPE \
        -g hs \
        -q 0.01 \
        --keep-dup all \
        --name "DUBR_${pool}" \
        --outdir "$peak_dir" \
        2> "$peak_dir/DUBR_${pool}_macs2.log"

    np="$peak_dir/DUBR_${pool}_peaks.narrowPeak"

    # Safety check: peaks and blacklist must both say "chr1" (or both say "1").
    # Otherwise bedtools would compare "1" with "chr1" and remove nothing.
    if [ "$(head -1 "$np" | cut -f1 | grep -c '^chr')" != "$(head -1 "$blacklist" | cut -f1 | grep -c '^chr')" ]; then
        echo "ERROR: chromosome names differ between peaks and blacklist - run 03b first" >&2
        exit 1
    fi

    # Column 7 of the peak file is the fold-enrichment. Keep rows where it's >= 5.
    awk -v fc="$fc_min" 'BEGIN{OFS="\t"} $7 >= fc' "$np" \
        > "$peak_dir/DUBR_${pool}_FC${fc_min}.narrowPeak"

    # Remove any peak that overlaps a blacklisted region (-v = keep non-overlaps).
    bedtools intersect -v \
        -a "$peak_dir/DUBR_${pool}_FC${fc_min}.narrowPeak" \
        -b "$blacklist" \
        > "$peak_dir/DUBR_${pool}_FC${fc_min}_noBL.narrowPeak"

    # Print how many peaks survived each filter (wc -l counts lines/peaks).
    n_fc=$(wc -l < "$peak_dir/DUBR_${pool}_FC${fc_min}.narrowPeak")
    n_bl=$(wc -l < "$peak_dir/DUBR_${pool}_FC${fc_min}_noBL.narrowPeak")
    echo "  $pool raw:        $(wc -l < "$np")"
    echo "  $pool FC>=$fc_min:     $n_fc"
    echo "  $pool FC+noBL:    $n_bl   (blacklist removed $((n_fc - n_bl)); expect a few dozen, not 0)"
done

echo "Peak files are in: $peak_dir"
