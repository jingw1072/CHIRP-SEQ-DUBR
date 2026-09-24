#!/usr/bin/env bash
# =============================================================================
# STEP 7 - Download the OTHER datasets the paper's figures are built on.
#
# Figures 5G-I, 6A, 6D-E and S8 overlay the DUBR peaks on:
#   * the authors' histone ChIP-seq (H3K4me3, H3K4me1, H3K27ac, H3K27me3 in
#     WT / dEX1 / dEX2 K562)          -> GEO GSE268500, ready-made bigWigs
#   * the authors' RNA-seq (WT/dEX1/dEX2) -> GEO GSE268502, bigWigs
#   * the authors' own ChIRP bigWigs (to compare with ours) -> GEO GSE268501
#   * public K562 TF ChIP-seq (HES1, GATA-1, NCOR1, RUNX1, TAF1) -> ENCODE
#   * UCSC hg19 gene annotation (refGene) for snapshots and TSS positions
# Total ~10 GB. Replicates are then averaged into one bigWig per mark x
# condition (what the paper plots). All files are hg19 with "chr" names.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 07_download_tracks.sh
# NOTE: needs outbound internet from the compute node. If downloads fail with
# a network error, run the script directly on the login node instead:
#     bash 07_download_tracks.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_tracks
#SBATCH --account=jiw619
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chirp-seq

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
trk="$out/tracks"
annot="$out/ref/annot"
mkdir -p "$trk"/{geo_chip,geo_rna,geo_chirp,encode,avg} "$annot"

geo_url() {   # geo_url GSM8292063 file.bigwig -> download URL
    local gsm=$1 file=$2
    echo "https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/${gsm}/suppl/${file}"
}
fetch() {     # fetch URL DEST  (skips files already downloaded)
    local url=$1 dest=$2
    if [ -s "$dest" ]; then echo "  have $(basename "$dest")"; return; fi
    echo "  get  $(basename "$dest")"
    curl -sSL --retry 5 --retry-delay 10 -o "$dest.part" "$url" && mv "$dest.part" "$dest"
}

# ---- 7a. Authors' histone ChIP-seq (GSE268500) -------------------------------
echo "=== GSE268500 histone ChIP-seq bigWigs ==="
chip=(
  GSM8292063_H3K4me3_WT_1  GSM8292064_H3K4me3_WT_2
  GSM8292065_H3K4me1_WT_1  GSM8292066_H3K4me1_WT_2
  GSM8292067_H3K27ac_WT_1  GSM8292068_H3K27ac_WT_2
  GSM8292069_H3K27me3_WT_1 GSM8292070_H3K27me3_WT_2
  GSM8292071_H3K4me3_EX1_1  GSM8292072_H3K4me3_EX1_2
  GSM8292073_H3K4me1_EX1_1  GSM8292074_H3K4me1_EX1_2
  GSM8292075_H3K27ac_EX1_1  GSM8292076_H3K27ac_EX1_2
  GSM8292077_H3K27me3_EX1_1 GSM8292078_H3K27me3_EX1_2
  GSM8628031_H3K4me3_EX2_1  GSM8628032_H3K4me3_EX2_2
  GSM8628033_H3K4me1_EX2_1  GSM8628034_H3K4me1_EX2_2
  GSM8628035_H3K27ac_EX2_1  GSM8628036_H3K27ac_EX2_2
  GSM8628037_H3K27me3_EX2_1 GSM8628038_H3K27me3_EX2_2
)
for f in "${chip[@]}"; do fetch "$(geo_url "${f%%_*}" "$f.bigwig")" "$trk/geo_chip/$f.bigwig"; done

# ---- 7b. Authors' RNA-seq (GSE268502) -----------------------------------------
echo "=== GSE268502 RNA-seq bigWigs ==="
rna=(
  GSM8292085_RNA_WT_1  GSM8292086_RNA_WT_2  GSM8292087_RNA_WT_3
  GSM8292088_RNA_EX1_1 GSM8292089_RNA_EX1_2 GSM8292090_RNA_EX1_3
  GSM8292091_RNA_EX2_1 GSM8292092_RNA_EX2_2 GSM8292093_RNA_EX2_3
)
for f in "${rna[@]}"; do fetch "$(geo_url "${f%%_*}" "$f.bigwig")" "$trk/geo_rna/$f.bigwig"; done

# ---- 7c. Authors' ChIRP bigWigs (GSE268501) -----------------------------------
echo "=== GSE268501 authors' ChIRP-seq bigWigs ==="
chirp=(
  GSM8292079_ChIRP_INPUT_1 GSM8292080_ChIRP_INPUT_2
  GSM8292081_ChIRP_EVEN_1  GSM8292082_ChIRP_EVEN_2
  GSM8292083_ChIRP_ODD_1   GSM8292084_ChIRP_ODD_2
)
for f in "${chirp[@]}"; do fetch "$(geo_url "${f%%_*}" "$f.bigwig")" "$trk/geo_chirp/$f.bigwig"; done

# ---- 7d. ENCODE K562 TF ChIP-seq (hg19, pooled fold-change over control) ------
echo "=== ENCODE K562 TF ChIP-seq ==="
enc() {   # enc ACCESSION LOCAL_NAME EXT   (ENCODE serves ".bigWig"; we store ".bw")
    local ext=$3 local_ext=$3
    [ "$ext" = "bigWig" ] && local_ext="bw"
    fetch "https://www.encodeproject.org/files/$1/@@download/$1.$ext" "$trk/encode/$2.$local_ext"
}
enc ENCFF071USV HES1_K562_FC   bigWig      # HES1  (GSE91470)
enc ENCFF833CBY GATA1_K562_FC  bigWig      # GATA1 (Snyder lab)
enc ENCFF211OZC NCOR1_K562_FC  bigWig      # NCOR1 (GSE92062)   Fig S8J
enc ENCFF962DJQ RUNX1_K562_FC  bigWig      # RUNX1 (GSE96253)   Fig S8J
enc ENCFF488POX TAF1_K562_FC   bigWig      # TAF1  (GSM803431)  Fig S8J
enc ENCFF676NPW HES1_K562_optimalIDR_peaks bed.gz   # HES1 peaks, Fig 6E

# ---- 7e. Average replicates -> one bigWig per mark x condition ---------------
# bigwigAverage is effectively single-threaded (~5-15 min per track), so six
# run at the same time in the background; `wait` pauses until they finish.
echo "=== averaging replicates (bigwigAverage, 6 at a time) ==="
avg() {   # avg OUT in1 in2 [in3]
    local o=$1; shift
    if [ -s "$trk/avg/$o.bw" ]; then echo "  have $o.bw"; return; fi
    echo "  make $o.bw"
    bigwigAverage -b "$@" -o "$trk/avg/$o.bw.part" -p 2 --binSize 50 \
        > "$trk/avg/$o.log" 2>&1 && mv "$trk/avg/$o.bw.part" "$trk/avg/$o.bw"
}
njobs=0
for mark in H3K4me3 H3K4me1 H3K27ac H3K27me3; do
    for cond in WT EX1 EX2; do
        r=( "$trk"/geo_chip/*_${mark}_${cond}_1.bigwig "$trk"/geo_chip/*_${mark}_${cond}_2.bigwig )
        avg "${mark}_${cond}" "${r[@]}" & njobs=$((njobs+1)); [ $((njobs % 6)) -eq 0 ] && wait
    done
done
for cond in WT EX1 EX2; do
    r=( "$trk"/geo_rna/*_RNA_${cond}_[123].bigwig )
    avg "RNA_${cond}" "${r[@]}" & njobs=$((njobs+1)); [ $((njobs % 6)) -eq 0 ] && wait
done
for pool in INPUT EVEN ODD; do
    r=( "$trk"/geo_chirp/*_ChIRP_${pool}_[12].bigwig )
    avg "ChIRPauthors_${pool}" "${r[@]}" & njobs=$((njobs+1)); [ $((njobs % 6)) -eq 0 ] && wait
done
wait
rm -f "$trk"/avg/*.part

# ---- 7f. Gene annotation: UCSC hg19 refGene ----------------------------------
if [ ! -s "$annot/refGene.txt.gz" ]; then
    echo "=== UCSC hg19 refGene ==="
    fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" "$annot/refGene.txt.gz"
fi
if [ ! -s "$annot/tss.bed" ]; then
    # one TSS per transcript: txStart on +, txEnd-1 on -  (BED: chr tss tss+1 gene . strand)
    zcat "$annot/refGene.txt.gz" \
      | awk 'BEGIN{OFS="\t"} $3 !~ /_/ { tss = ($4=="+") ? $5 : $6-1; print $3, tss, tss+1, $13, ".", $4 }' \
      | sort -k1,1 -k2,2n -u > "$annot/tss.bed"
    echo "  tss.bed: $(wc -l < "$annot/tss.bed") transcript TSSs"
fi

# ---- 7g. Sanity check: every track must be hg19 (chr1 = 249,250,621 bp) ------
echo "=== assembly check ==="
python - "$trk" <<'EOF'
import sys, glob, pyBigWig
trk = sys.argv[1]; bad = 0
for f in sorted(glob.glob(f"{trk}/avg/*.bw") + glob.glob(f"{trk}/encode/*.bw")):
    bw = pyBigWig.open(f); l = bw.chroms().get("chr1"); bw.close()
    tag = "hg19" if l == 249250621 else ("hg38?" if l == 248956422 else f"chr1={l}")
    bad += tag != "hg19"
    print(f"  {tag:6s} {f.split('/')[-1]}")
sys.exit(1 if bad else 0)
EOF
echo "Tracks are in: $trk   (averaged: $trk/avg)"
