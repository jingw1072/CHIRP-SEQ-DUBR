#!/usr/bin/env bash
#SBATCH --account=htl145           # our allocation (who gets billed)
#SBATCH --partition=hotel          # which set of machines to run on
#SBATCH --qos=hotel                # quality-of-service that matches the partition
#SBATCH --nodes=1                  # use 1 computer
#SBATCH --ntasks=1                 # run 1 task
#SBATCH --cpus-per-task=8          # give that task 8 CPU cores
#SBATCH --mem=16G                  # give it 16 GB of memory
#SBATCH --time=12:00:00            # kill the job if it runs longer than 12 hours

# ---- Turn on our software (conda environment with all the tools) ------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Where everything lives -------------------------------------------------
# Big files go on Lustre scratch (the fast, large disk), NOT your home folder.
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
fastq_dir="$out/fastq"      # final renamed FASTQ files end up here
sra_dir="$out/sra"          # raw downloaded .sra files
tmp_dir="$out/tmp"          # scratch space used while converting sra -> fastq
ref_dir="$out/ref"          # genome index + blacklist
threads=8

# Make the folders (the -p means "don't complain if they already exist").
mkdir -p "$fastq_dir" "$sra_dir" "$tmp_dir" "$ref_dir"

# ---- Download the 6 samples -------------------------------------------------
# Each line is: <SRR accession><TAB><friendly name>. The friendly name tells us
# the role (Input / EVEN / ODD) and which replicate (rep1 / rep2).
samples=(
  "SRR29207792	Input_rep1"
  "SRR29207791	Input_rep2"
  "SRR29207790	EVEN_rep1"
  "SRR29207789	EVEN_rep2"
  "SRR29207788	ODD_rep1"
  "SRR29207787	ODD_rep2"
)

for entry in "${samples[@]}"; do
    # Split the line into the SRR id and the friendly name.
    srr="${entry%%$'\t'*}"     # text before the tab
    name="${entry##*$'\t'}"    # text after the tab
    echo "=== ${name} (${srr}) ==="

    # 1) Download the raw .sra file from NCBI.
    prefetch --max-size 100G -O "$sra_dir" "$srr"

    # 2) Convert .sra to FASTQ. --split-files makes _1 and _2 (paired-end reads).
    fasterq-dump --split-files --threads "$threads" \
        --temp "$tmp_dir" \
        --outdir "$fastq_dir" "$sra_dir/$srr/$srr.sra"

    # 3) Compress and rename to our friendly names (R1 = read 1, R2 = read 2).
    gzip -f "$fastq_dir/${srr}_1.fastq"
    gzip -f "$fastq_dir/${srr}_2.fastq"
    mv -f "$fastq_dir/${srr}_1.fastq.gz" "$fastq_dir/${name}_R1.fastq.gz"
    mv -f "$fastq_dir/${srr}_2.fastq.gz" "$fastq_dir/${name}_R2.fastq.gz"
    echo "Done: ${name}"
done

echo "All 6 FASTQ files are in: $fastq_dir"

# ---- Download the reference genome (hg19) and the blacklist -----------------
# The genome index lets bowtie2 align reads. We only download it if it's missing.
if [ ! -e "$ref_dir/GRCh37.1.bt2" ] && [ ! -e "$ref_dir/GRCh37.1.bt2l" ]; then
    echo "=== downloading hg19 bowtie2 index ==="
    cd "$ref_dir"
    wget https://genome-idx.s3.amazonaws.com/bt/GRCh37.zip
    unzip -o GRCh37.zip
    rm -f GRCh37.zip
    cd -
fi

# The blacklist marks problem regions of the genome that we throw peaks out of.
if [ ! -e "$ref_dir/hg19-blacklist.v2.bed" ]; then
    echo "=== downloading hg19 ENCODE blacklist ==="
    wget -O "$ref_dir/hg19-blacklist.v2.bed.gz" \
        https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
    gunzip -f "$ref_dir/hg19-blacklist.v2.bed.gz"
fi

echo "Reference genome is ready in: $ref_dir"
