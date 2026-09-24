#!/usr/bin/env bash
# =============================================================================
# STEP 12 - Motif analysis of the DUBR ChIRP-seq peaks.
#
# WHAT THE PAPER DID
#   * HOMER findMotifsGenome.pl on *ChIP-seq* peaks (Fig 4C: gained / reduced
#     H3K27ac; Fig S8I: HES1 peaks). NOT on the ChIRP peaks.
#   * On the ChIRP peaks, one motif step that is only visible in their R code
#     ("ChIRP peaks vs DREME"): DREME on the probe sequences, then peaks that
#     contain a probe-like motif were flagged as suspect and the remainder kept
#     as "TRUE peaks" for the heatmaps. That is a ChIRP-specific QC step.
#   * Fig 5F / 6F / S8A are gene-ontology enrichment (Metascape), not motifs
#     -> step 13.
#
# WHAT A STANDARD ChIRP-seq MOTIF ANALYSIS LOOKS LIKE (Chu et al. 2011 and
# most lncRNA ChIRP papers since) - all three parts are done here:
#   A. Probe off-target check (ChIRP-specific). Biotinylated probes can
#      hybridise to genomic DNA directly, giving peaks that have nothing to do
#      with the RNA. Two ways to find them:
#        A1 align the probe sequences to the genome  -> peaks overlapping a hit
#        A2 look for probe FRAGMENTS (exact k-mers) inside peak sequences - the
#           authors' DREME-based idea, made explicit
#      Expected: the DUBR gene itself (probes bind the nascent RNA there AND
#      match the DNA), plus a handful of off-targets.
#   B. De novo motif discovery in the peaks (HOMER). Does the RNA occupy a
#      particular DNA sequence? (HOTAIR: GA-rich polypurine motif; many
#      lncRNAs: none, or motifs of the TF that recruits them.)
#   C. Known TF-motif enrichment in the peaks (HOMER known motifs vs a
#      GC-matched genomic background). Tells you which TFs co-occupy DUBR
#      sites - candidate protein partners / recruiters.
#   (D. RNA:DNA triplex prediction (Triplexator, TDF) is another common lncRNA
#      analysis; it needs the DUBR transcript sequence and is not covered.)
#
# INPUT  peaks/DUBR_peaks_final.bed (step 06), peaks/DUBR_EVEN_peaks.narrowPeak
#        (step 04, for the summit position), ref/hg19.fa, bowtie2 index (step 00)
# OUTPUT motifs/  (see the echo lines at the end)
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 12_motifs_chirp.sh
# Needs the chirp-motif environment (README section 2c). Downloads hg19.fa from
# UCSC on first use (~900 MB) - needs internet on the node. Three HOMER runs:
# allow ~2 h.
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=dubr_motifs
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/dubr/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chirp-motif        # homer, bedtools, samtools, bowtie2

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr"
D="$out"
peak_dir="$D/peaks"
mot="$D/motifs"
genome_fa="$D/ref/hg19.fa"                   # UCSC hg19, "chr" names (downloaded below)
bt2_index="$D/ref/GRCh37"                    # bowtie2 index from step 00 (Ensembl names)
probes="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}/dubr_chirp_probes.fa"   # Table S6, next to this script

# hg19 FASTA for sequence extraction and HOMER (once)
if [ ! -s "$genome_fa" ]; then
    echo "=== downloading UCSC hg19.fa (~900 MB) ==="
    curl -sSL --retry 5 -o "$genome_fa.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gunzip -f "$genome_fa.gz" && samtools faidx "$genome_fa"
fi
threads=16
mkdir -p "$mot"

# ---- 12a. Sequences to analyse ----------------------------------------------
# Motif tools want short, centred windows: +/-100 bp around the MACS2 summit
# (narrowPeak column 10 = summit offset from the peak start). Our final peaks
# are EVEN intervals, so their summits are in the EVEN narrowPeak file.
echo "=== 12a. summit-centred 200 bp windows for $(wc -l < "$peak_dir/DUBR_peaks_final.bed") final peaks ==="
awk 'BEGIN{OFS="\t"}
     NR==FNR { summit[$4] = $2 + $10; next }          # 1st file: name -> summit position
     ($4 in summit) { s = summit[$4]; print $1, s-100, s+100, $4 }' \
    "$peak_dir/DUBR_EVEN_peaks.narrowPeak" "$peak_dir/DUBR_peaks_final.bed" \
    | sort -k1,1 -k2,2n > "$mot/peaks_summit200.bed"
bedtools getfasta -fi "$genome_fa" -bed "$peak_dir/DUBR_peaks_final.bed" -nameOnly -fo "$mot/peaks_full.fa"
bedtools getfasta -fi "$genome_fa" -bed "$mot/peaks_summit200.bed"      -nameOnly -fo "$mot/peaks_summit200.fa"

# ---- 12b. A1: align the probes to the genome ---------------------------------
# 20-mers; report up to 500 alignments per probe with <=1 mismatch (bowtie2
# mismatch penalty 6 -> --score-min L,-6,0). Seeds of 10 bp at every offset
# (-L 10 -i S,1,0) guarantee every <=1-mismatch hit is found; allowing seed
# mismatches (-N 1) or exhaustive -a makes this take forever on a human genome.
# Looser matches inside peaks are covered by the FIMO scan in 12c. The index
# uses Ensembl names ("1"), so add "chr" when converting to BED.
echo "=== 12b. probe off-target check: A1 bowtie2 alignment of probes (<=1 mismatch) ==="
bowtie2 -f -k 500 --end-to-end -L 10 -N 0 -i S,1,0 --score-min L,-6,0 \
    -x "$bt2_index" -U "$probes" -p 4 --no-unal 2> "$mot/probe_bowtie2.log" \
  | samtools view - \
  | awk 'BEGIN{OFS="\t"} { chrom = ($3 ~ /^chr/) ? $3 : "chr"$3; if (chrom=="chrMT") chrom="chrM"
                           print chrom, $4-1, $4-1+length($10), $1, ($2 % 32 >= 16) ? "-" : "+" }' \
  | sort -k1,1 -k2,2n > "$mot/probe_genome_hits.bed"
echo "  probe alignments (<=1 mismatch): $(wc -l < "$mot/probe_genome_hits.bed")"
cut -f4 "$mot/probe_genome_hits.bed" | sort | uniq -c | awk '{printf "    %-8s %d hits\n", $2, $1}'
bedtools intersect -u -a "$peak_dir/DUBR_peaks_final.bed" -b "$mot/probe_genome_hits.bed" \
    > "$mot/peaks_with_probe_alignment.bed"
echo "  final peaks overlapping a probe alignment: $(wc -l < "$mot/peaks_with_probe_alignment.bed")"
cat "$mot/peaks_with_probe_alignment.bed" | sed 's/^/    /'

# ---- 12c. A2: probe k-mer scan inside the peaks --------------------------------
# A full-length alignment (A1) only catches near-perfect 20-nt matches. But a
# probe can also hybridise through a PARTIAL match, and de novo motif discovery
# on these peaks (12d, below) rediscovers the probe sequences themselves as the
# top motifs - i.e. many peaks contain probe fragments. This is the effect the
# authors' DREME step was after. First, the transparent version: does a peak
# contain any EXACT k-mer of a probe (either strand)?  Pure awk + grep:
#   1. all k-mers of the 6 DUBR probes, plus reverse complements -> probe_kmers_K.txt
#   2. peaks_full.fa linearised to "name<TAB>sequence"
#   3. grep -F -f (fixed strings from file)  -> names of peaks with a hit
# (shuffled-sequence control: k=10 5%, k=11 0.6%, k>=12 0 - so every exact
#  12-mer hit is real).
echo "=== 12c. probe off-target check: A2 exact probe k-mer scan of peak sequences ==="
grep -v '^>' "$probes" | grep -v '^$' | head -6 > "$mot/probe_seqs.txt"          # the 6 DUBR probes
awk '/^>/ { if (n) print n "\t" seq; n = substr($1, 2); seq = ""; next }
     { seq = seq toupper($0) }  END { print n "\t" seq }' "$mot/peaks_full.fa" > "$mot/peaks_full.tsv"
kmers() {   # kmers K  -> all K-mers of the probes, both strands, one per line
    awk -v k="$1" 'function rc(x,  i, r, c) { r = ""; for (i = length(x); i >= 1; i--) { c = substr(x, i, 1);
                       r = r (c=="A"?"T":c=="T"?"A":c=="C"?"G":c=="G"?"C":c) } return r }
         { for (i = 1; i <= length($1)-k+1; i++) { m = substr($1, i, k); print m; print rc(m) } }' "$mot/probe_seqs.txt" | sort -u
}
echo "  k    peaks containing an exact probe k-mer"
for k in 10 11 12 13 14 16 18 20; do
    kmers "$k" > "$mot/probe_kmers_$k.txt"
    n=$(grep -F -f "$mot/probe_kmers_$k.txt" "$mot/peaks_full.tsv" | wc -l)
    awk -v k="$k" -v n="$n" -v t="$(wc -l < "$mot/peaks_full.tsv")" 'BEGIN{printf "  %-4s %5d  (%4.1f%%)\n", k, n, 100*n/t}'
done
# Exact matching under-counts: HOMER de novo discovery (12d) finds the probes
# even in the exact-12-mer-free peaks, with one mismatch. So the flag uses a
# mismatch-TOLERANT scan, HOMER-native: seq2profile.pl turns every probe 14-mer
# into a motif allowing 1 mismatch; homer2 find scans both strands of every peak.
# Calibration (this data; shuffled-sequence control in brackets):
#   12-mer <=1 mm: 88% (8.3%)   14-mer <=1 mm: 77% (0.2%)   16-mer <=1 mm: 42% (0)
#   14-mer <=2 mm: 89% (8.3%)   16-mer <=2 mm: 71% (0.4%)
# -> 14-mer with <=1 mismatch: many flagged, essentially no chance hits.
K=14; MM=1
> "$mot/probe_${K}mer_${MM}mm.motifs"; i=0
while read -r p; do
    for ((s = 0; s <= ${#p} - K; s++)); do
        i=$((i + 1)); seq2profile.pl "${p:$s:$K}" "$MM" "probe_frag_$i" >> "$mot/probe_${K}mer_${MM}mm.motifs" 2>/dev/null
    done
done < "$mot/probe_seqs.txt"
homer2 find -i "$mot/peaks_full.fa" -m "$mot/probe_${K}mer_${MM}mm.motifs" -p "$threads" 2>/dev/null \
    | cut -f1 | sort -u > "$mot/peaks_with_probe_${K}mer_${MM}mm.txt"
echo "  -> any probe ${K}-mer with <=${MM} mismatch (either strand): $(wc -l < "$mot/peaks_with_probe_${K}mer_${MM}mm.txt") peaks"
# Flag = exact 12-mer OR 14-mer with <=1 mismatch (both criteria have ~0 chance
# hits but catch different fragments: a 12-mer with mismatched flanks fails the
# 14-mer test). Residual shorter/degenerate fragments still exist - see the de
# novo motifs of the clean set - so "clean" means "no strong probe match", not
# "no probe sequence at all".
grep -F -f "$mot/probe_kmers_12.txt" "$mot/peaks_full.tsv" | cut -f1 \
    | cat - "$mot/peaks_with_probe_${K}mer_${MM}mm.txt" | sort -u > "$mot/peaks_with_probe_kmer.txt"
echo "  -> flagged (exact 12-mer OR ${K}-mer <=${MM} mm): $(wc -l < "$mot/peaks_with_probe_kmer.txt") peaks"

# "TRUE peaks" (authors' term) = peaks with neither a probe alignment nor a probe fragment
cat "$mot/peaks_with_probe_kmer.txt" <(cut -f4 "$mot/peaks_with_probe_alignment.bed") | sort -u > "$mot/peaks_probe_flagged.txt"
awk 'NR==FNR{bad[$1]=1; next} !($4 in bad)' "$mot/peaks_probe_flagged.txt" "$peak_dir/DUBR_peaks_final.bed" \
    > "$peak_dir/DUBR_peaks_final_noProbe.bed"
echo "  flagged by A1 or A2: $(wc -l < "$mot/peaks_probe_flagged.txt")  ->  clean set: $(wc -l < "$peak_dir/DUBR_peaks_final_noProbe.bed") peaks (peaks/DUBR_peaks_final_noProbe.bed)"
# Are the biologically interesting peaks affected? (C1 cluster from step 08; peaks of the direct targets from step 09)
for f in "$peak_dir/clusters/C1.bed" "$peak_dir/DUBR_peaks_direct_upregulated.bed" "$peak_dir/DUBR_peaks_direct_downregulated.bed"; do
    [ -s "$f" ] || continue
    n_all=$(wc -l < "$f")
    n_clean=$(bedtools intersect -u -f 1.0 -r -a "$f" -b "$peak_dir/DUBR_peaks_final_noProbe.bed" | wc -l)
    printf "  %-40s %4d peaks, %4d without probe sequence\n" "$(basename "$f")" "$n_all" "$n_clean"
done

# ---- 12d. B + C: HOMER de novo discovery and known-motif enrichment -----------
# findMotifsGenome.pl picks GC-matched random genomic regions as background,
# then (i) tests known vertebrate TF motifs (knownResults.txt) and (ii) searches
# for new motifs of length 8/10/12 (homerResults.html).
# -size given = use our 200 bp windows as they are. -preparsedDir keeps HOMER's
# genome cache out of the reference folder. -mset vertebrates restricts the
# known-motif library to vertebrate TFs (the default mixes in plant/yeast/fly
# motifs). -S 15 = report up to 15 de novo motifs. Each run takes ~30-45 min.
#
# Run 1: ALL final peaks. Expect the de novo motifs to be probe fragments
#        (the point of 12c). Run 2: the clean set - the TF analysis that matters.
homer_run() {   # homer_run BED OUTDIR [extra args]
    local bed=$1 outdir=$2; shift 2
    if [ -s "$outdir/knownResults.txt" ]; then echo "  (already done - delete $outdir to redo)"; return; fi
    findMotifsGenome.pl "$bed" "$genome_fa" "$outdir" -size given -mset vertebrates -S 15 \
        -p "$threads" -preparsedDir "$mot/homer_preparsed" "$@" > "$outdir.log" 2>&1
}
show_known() { awk -F'\t' -v n="$2" 'NR>1 && NR<=n+1 {printf "    %-45s p=%-8s %s of peaks (bg %s)\n", substr($1,1,45), $3, $7, $9}' "$1/knownResults.txt"; }
show_denovo() {   # first 8 de novo motifs (no "| head": under pipefail an early-closing head aborts the script)
    for i in 1 2 3 4 5 6 7 8; do
        f="$1/homerResults/motif$i.motif"; [ -f "$f" ] || continue
        head -1 "$f" | awk -F'\t' '{ sub(/^>/,"",$1); split($2,a,"BestGuess:"); split($6,b,",")
                                    printf "    %-14s best match %-40s %s (bg %s)\n", $1, substr(a[2],1,40), b[1], b[2] }'
    done
}
summits() {   # summits IN.bed OUT.bed : +/-100 bp around the MACS2 summit; names recovered by exact coordinates
    bedtools intersect -u -f 1.0 -r -a "$peak_dir/DUBR_peaks_final.bed" -b "$1" \
      | awk 'BEGIN{OFS="\t"} NR==FNR{summit[$4]=$2+$10; next} ($4 in summit){s=summit[$4]; print $1,s-100,s+100,$4}' \
            "$peak_dir/DUBR_EVEN_peaks.narrowPeak" - | sort -k1,1 -k2,2n > "$2"
}

echo "=== 12d. HOMER, run 1: all $(wc -l < "$mot/peaks_summit200.bed") final peaks ==="
homer_run "$mot/peaks_summit200.bed" "$mot/homer_final_peaks"
echo "  de novo motifs (compare with the probes in $probes):"; show_denovo "$mot/homer_final_peaks"
echo "  top known motifs:"; show_known "$mot/homer_final_peaks" 8

summits "$peak_dir/DUBR_peaks_final_noProbe.bed" "$mot/noProbe_summit200.bed"
echo "=== 12d. HOMER, run 2: $(wc -l < "$mot/noProbe_summit200.bed") peaks without probe sequence ==="
homer_run "$mot/noProbe_summit200.bed" "$mot/homer_noProbe_peaks"
echo "  de novo motifs:"; show_denovo "$mot/homer_noProbe_peaks"
echo "  top known motifs:"; show_known "$mot/homer_noProbe_peaks" 12

# ---- 12e. Known motifs in the regulatory-element cluster C1 (step 08) --------
# Small set (~70 peaks): known-motif enrichment only (-nomotif), suggestive at best.
if [ -s "$peak_dir/clusters/C1.bed" ]; then
    summits "$peak_dir/clusters/C1.bed" "$mot/C1_summit200.bed"
    echo "=== 12e. HOMER known motifs, cluster C1 ($(wc -l < "$mot/C1_summit200.bed") peaks) ==="
    homer_run "$mot/C1_summit200.bed" "$mot/homer_C1" -nomotif
    show_known "$mot/homer_C1" 8
fi

echo
echo "outputs in $mot :"
echo "  probe_genome_hits.bed, peaks_with_probe_alignment.bed  (A1)"
echo "  probe_kmers_*.txt (exact sweep), probe_14mer_1mm.motifs, peaks_with_probe_kmer.txt  (A2)"
echo "  homer_final_peaks/    all peaks   : homerResults.html (de novo = probe fragments), knownResults.html"
echo "  homer_noProbe_peaks/  clean peaks : homerResults.html, knownResults.html  <- the TF analysis to read"
echo "  homer_C1/             C1 cluster  : knownResults.html"
echo "  ../peaks/DUBR_peaks_final_noProbe.bed  = the authors' \"TRUE peaks\" equivalent"
