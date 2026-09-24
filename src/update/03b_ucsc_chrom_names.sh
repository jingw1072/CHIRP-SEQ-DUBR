#!/usr/bin/env bash
# =============================================================================
# STEP 3b - Rename chromosomes in the merged BAMs to UCSC style (chr1, chrX...).
#
# WHY: the genome index from step 00 (GRCh37.zip) names chromosomes "1", "2",
# "X"... but the blacklist (step 04), the authors' GEO tracks, ENCODE tracks and
# GREAT all use "chr1", "chr2", "chrX". `bedtools intersect` treats "1" and
# "chr1" as DIFFERENT chromosomes, so without this step the blacklist filter in
# step 04 silently removes nothing. Only the BAM *header* is rewritten (fast).
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 03b_ucsc_chrom_names.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_chrnames
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/dubr/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
bam_dir="$out/bam"
threads=8

for cond in Input EVEN ODD; do
    in="$bam_dir/${cond}_merged.bam"
    outbam="$bam_dir/${cond}_merged.chr.bam"

    # If the index already used "chr" names there is nothing to do: just link.
    if samtools view -H "$in" | grep -q $'^@SQ\tSN:chr1\t'; then
        echo "$cond: already UCSC-named; linking"
        ln -sf "$in" "$outbam"; ln -sf "$in.bai" "$outbam.bai"
        continue
    fi

    echo "=== reheader $cond -> UCSC names ==="
    # 1 -> chr1 ... X -> chrX, Y -> chrY, MT -> chrM (other contigs untouched)
    samtools view -H "$in" \
      | sed -E 's/^(@SQ\tSN:)([0-9]+|X|Y)\t/\1chr\2\t/; s/^(@SQ\tSN:)MT\t/\1chrM\t/' \
      > "$bam_dir/${cond}_merged.chr.header.sam"
    samtools reheader "$bam_dir/${cond}_merged.chr.header.sam" "$in" > "$outbam"
    samtools index -@ "$threads" "$outbam"
    rm -f "$bam_dir/${cond}_merged.chr.header.sam"
    echo "  $(samtools view -H "$outbam" | grep -c '^@SQ.SN:chr') chr-named contigs"
done

echo "UCSC-named BAMs are in: $bam_dir (*_merged.chr.bam)"
