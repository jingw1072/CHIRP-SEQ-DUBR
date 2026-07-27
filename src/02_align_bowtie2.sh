#!/usr/bin/env bash
# =============================================================================
# STEP 2 - Align the cleaned reads to the human genome (hg19) with bowtie2.
#
# "Alignment" means figuring out where in the genome each short read came from.
# The output is a BAM file (a compressed table of aligned reads), sorted by
# genome position so later tools can use it quickly.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 02_align_bowtie2.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_align
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24         # hotel nodes allow up to 28 cores per node
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
fastq_dir="$out/fastq_trimmed"    # input: cleaned reads from step 01b
bam_dir="$out/bam"                # output: aligned BAM files
index="$out/ref/GRCh37"           # the genome index downloaded in step 00
bt2_threads=16                    # cores for bowtie2 (the aligning)
sort_threads=6                    # cores for samtools sort (the sorting)
sort_mem="4G"                     # memory per sort thread (6 x 4G = 24G)

mkdir -p "$bam_dir"

for sample in Input_rep1 Input_rep2 EVEN_rep1 EVEN_rep2 ODD_rep1 ODD_rep2; do
    echo "=== aligning $sample ==="
    # bowtie2 aligns the reads and prints results; we pipe ("|") those results
    # straight into "samtools sort" so we never write the huge unsorted file.
    #   -x   the genome index
    #   -1/-2  the two paired read files
    #   -X 1000  allow read pairs up to 1000 bp apart
    #   --no-unal  drop reads that didn't align anywhere
    #   --mm  share the genome index in memory across samples (faster)
    bowtie2 \
        -x "$index" \
        -1 "$fastq_dir/${sample}_R1.fastq.gz" \
        -2 "$fastq_dir/${sample}_R2.fastq.gz" \
        -X 1000 \
        --no-unal \
        --mm \
        -p "$bt2_threads" \
        2> "$bam_dir/${sample}_bowtie2.log" \
      | samtools sort -@ "$sort_threads" -m "$sort_mem" \
            -T "$bam_dir/${sample}_sorttmp" \
            -o "$bam_dir/${sample}_sorted.bam" -

    # Build an index for the BAM (a .bai file) so tools can jump around it fast.
    samtools index -@ "$sort_threads" "$bam_dir/${sample}_sorted.bam"
    echo "Done: $sample"
done

echo "Aligned BAM files are in: $bam_dir"
