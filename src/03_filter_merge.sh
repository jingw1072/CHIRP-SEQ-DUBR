#!/usr/bin/env bash
# =============================================================================
# STEP 3 - Clean up the aligned reads, then combine the two replicates.
#
# Two things happen here:
#   1) FILTER each sample: remove PCR duplicates, keep only properly-paired
#      reads that aligned well and uniquely.
#   2) MERGE: combine rep1 + rep2 of each condition (Input / EVEN / ODD) into
#      one BAM, because the next step calls peaks on the combined data.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 03_filter_merge.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_filter
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
bam_dir="$out/bam"
threads=16
sort_mem="2G"
tmp="$bam_dir/sorttmp"   # temporary files for sorting live on scratch

# ---- 1) Filter each sample one at a time ------------------------------------
for sample in Input_rep1 Input_rep2 EVEN_rep1 EVEN_rep2 ODD_rep1 ODD_rep2; do
    echo "=== filtering $sample ==="
    # To find PCR duplicates in paired-end data, samtools needs the reads in a
    # specific order and with mate info attached. So we go through these steps:
    #   sort by name -> fixmate (add mate info) -> sort by position -> markdup
    samtools sort -n -@ "$threads" -m "$sort_mem" -T "${tmp}_${sample}_n" "$bam_dir/${sample}_sorted.bam" -o "$bam_dir/${sample}_namesort.bam"
    samtools fixmate -m -@ "$threads" "$bam_dir/${sample}_namesort.bam" "$bam_dir/${sample}_fixmate.bam"
    samtools sort -@ "$threads" -m "$sort_mem" -T "${tmp}_${sample}_p" "$bam_dir/${sample}_fixmate.bam" -o "$bam_dir/${sample}_possort.bam"
    # -r actually removes the duplicates it finds.
    samtools markdup -r -@ "$threads" "$bam_dir/${sample}_possort.bam" "$bam_dir/${sample}_dedup.bam" \
        2> "$bam_dir/${sample}_markdup.log"

    # Now keep only the good reads:
    #   -f 2     keep properly-paired reads
    #   -F 1280  drop duplicates + secondary alignments
    #   -q 30    drop poorly/low-quality mapped reads
    samtools view -b -@ "$threads" -f 2 -F 1280 -q 30 \
        "$bam_dir/${sample}_dedup.bam" > "$bam_dir/${sample}_filt.bam"
    samtools index -@ "$threads" "$bam_dir/${sample}_filt.bam"

    # Delete the intermediate files we no longer need (saves disk space).
    rm -f "$bam_dir/${sample}_namesort.bam" "$bam_dir/${sample}_fixmate.bam" "$bam_dir/${sample}_possort.bam" "$bam_dir/${sample}_dedup.bam"
done

# ---- 2) Merge the two replicates of each condition --------------------------
for cond in Input EVEN ODD; do
    echo "=== merging $cond replicates ==="
    samtools merge -f -@ "$threads" "$bam_dir/${cond}_merged.bam" \
        "$bam_dir/${cond}_rep1_filt.bam" "$bam_dir/${cond}_rep2_filt.bam"
    samtools index -@ "$threads" "$bam_dir/${cond}_merged.bam"
    # flagstat writes a small summary (how many reads) to a text file.
    samtools flagstat "$bam_dir/${cond}_merged.bam" > "$bam_dir/${cond}_merged_flagstat.txt"
done

echo "Filtered and merged BAM files are in: $bam_dir"
