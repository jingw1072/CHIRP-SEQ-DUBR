#!/usr/bin/env bash
# =============================================================================
# STEP 10 - Genome-browser snapshots with pyGenomeTracks:
#           Fig 5I (ATF5, TAF12, CSF1, ACVRL1), Fig 6A (HES1), DUBR-locus check.
#
# HOW pyGenomeTracks WORKS (the whole idea in two lines)
#   1. a plain-text .ini file lists the tracks, top to bottom   -> tracks_ini/*.ini
#   2. one command draws any region:  pyGenomeTracks --tracks X.ini --region chr:a-b -o out.png
#
# The .ini files live next to this script in tracks_ini/ and are meant to be
# read and edited by hand (colour, height, max_value ...). Their file paths are
# relative to the analysis folder, so this script `cd`s there before drawing.
# `make_tracks_file --trackFiles a.bw b.bed -o new.ini` writes a template .ini
# if you want to start a new figure.
#
# The gene track is UCSC refGene (downloaded in step 07) converted to BED12,
# one (longest) transcript per gene so the panel stays readable.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 10_tracks_pygenometracks.sh            # all figures
#     sbatch 10_tracks_pygenometracks.sh fig6A      # one figure
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_tracks_fig
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

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chirp-fig                     # pyGenomeTracks (README section 2b)

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
D="$out"                                     # analysis folder (bigwig/, tracks/, peaks/)
refgene="$D/ref/annot/refGene.txt.gz"        # from step 07
ini_dir="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}/tracks_ini"   # folder you ran sbatch from
fig_dir="$D/figures/tracks"
which=${1:-all}
mkdir -p "$fig_dir"

# ---- gene models: refGene -> BED12, longest transcript per gene (once) ------
# refGene columns: 3 chrom, 4 strand, 5 txStart, 6 txEnd, 7 cdsStart, 8 cdsEnd,
# 9 exonCount, 10 exonStarts, 11 exonEnds, 13 gene name.
# BED12 = chrom start end name score strand thickStart thickEnd rgb nExons sizes starts
genes="$D/tracks/refGene_hg19_genes.bed"
if [ ! -s "$genes" ]; then
    echo "=== building gene BED12 from refGene ==="
    zcat "$refgene" \
      | awk 'BEGIN{OFS="\t"} $3 !~ /_/ {
              n = split($10, s, ","); split($11, e, ",")
              sizes = ""; starts = ""
              for (i = 1; i < n; i++) { sizes = sizes (e[i]-s[i]) ","; starts = starts (s[i]-$5) "," }
              print $3, $5, $6, $13, 0, $4, $7, $8, 0, $9, sizes, starts, $6-$5 }' \
      | sort -k4,4 -k13,13nr | awk '!seen[$4]++' \
      | cut -f1-12 | sort -k1,1 -k2,2n > "$genes"
    echo "  $(wc -l < "$genes") genes"
fi

# ---- figures: <ini name>  <locus label>  <hg19 region> ----------------------
# Regions match the paper's snapshots (Fig 5I, Fig 6A) + the DUBR gene itself.
snapshots=(
  "fig5I       ATF5        chr19:50393938-50437327"
  "fig5I       TAF12       chr1:28926824-28990457"
  "fig5I       CSF1        chr1:110412575-110486118"
  "fig5I       ACVRL1      chr12:52287848-52321459"
  "fig6A       HES1        chr3:193832080-193860255"
  "dubr_locus  DUBR_locus  chr3:106800000-107070000"
)

cd "$D"                                      # .ini paths are relative to here
for entry in "${snapshots[@]}"; do
    read -r ini label region <<< "$entry"
    [ "$which" = "all" ] || [ "$which" = "$ini" ] || continue
    echo "=== $ini  $label  $region ==="
    pyGenomeTracks \
        --tracks "$ini_dir/$ini.ini" \
        --region "$region" \
        --outFileName "$fig_dir/${ini}_${label}.png" \
        --width 18 --dpi 200 --fontSize 9 --trackLabelFraction 0.22 \
        2> >(grep -iE "error|warn" >&2 || true)
done
echo "snapshots in $fig_dir"
