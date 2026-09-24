# DUBR ChIRP-seq analysis on the UCSD TSCC cluster

This folder runs the DUBR ChIRP-seq analysis (https://pmc.ncbi.nlm.nih.gov/articles/PMC11850227) on UC San Diego's supercomputer,
the **Triton Shared Computing Cluster (TSCC)**. You submit the scripts **one at
a time**, in numbered order (00, 01, 01b, 02, 03, 04, 05, 06). You run it, wait for it to finish, then run the next.

The scripts are set up for the `jiw619` allocation on the `hotel`
partition.

---

## 1. Log in to TSCC

From your own computer's terminal:

```bash
ssh <your_username>@login.tscc.sdsc.edu
```

Type your UCSD AD password, then approve the DUO 2-step prompt on your phone.

---

## 2. One-time setup: create the software environment

Now create a conda environment named
`chirp-seq` that contains every tool the scripts use. Run these commands **once**
on a login node:

```bash

# Create the environment with all the tools (this can take a few minutes)
conda create -y -n chirp-seq -c bioconda -c conda-forge \
    sra-tools \
    fastqc \
    multiqc \
    fastp \
    bowtie2 \
    samtools \
    bedtools \
    macs2 \
    deeptools

# Check it works
conda activate chirp-seq
which bowtie2 macs2   # should print paths inside .../chirp-seq/bin/
```

The scripts turn this environment on automatically (`conda activate chirp-seq`),
so you don't have to activate it yourself before submitting.

### 2b. Second environment for the figure scripts (steps 09 and 10)

Steps 09 (R: GREAT, Venns, heatmap) and 10 (pyGenomeTracks snapshots) use a
separate environment called `chirp-fig`, so the R/Bioconductor packages don't
disturb the analysis tools above. Create it once on a login node (10–20 min;
`mamba` is faster than `conda` if you have it):

```bash
# conda's package cache can be large; keep it on scratch, not in your 100 GB home
export CONDA_PKGS_DIRS=/tscc/lustre/ddn/scratch/$USER/conda_pkgs

conda create -y -n chirp-fig -c conda-forge -c bioconda \
    'r-base>=4.3,<4.5' \
    bioconductor-rgreat bioconductor-txdb.hsapiens.ucsc.hg19.knowngene bioconductor-org.hs.eg.db \
    bioconductor-complexheatmap bioconductor-rtracklayer bioconductor-genomicranges \
    r-ggplot2 r-eulerr r-data.table r-ggrepel r-dplyr r-tidyr r-ggpubr r-circlize r-readxl \
    pygenometracks

# Check it works
conda activate chirp-fig
Rscript -e 'library(rGREAT); library(ComplexHeatmap)' && pyGenomeTracks --version
```

### 2c. Third environment for motif analysis (step 12)

HOMER has many dependencies of its own, so it gets its own environment:

```bash
conda create -y -n chirp-motif -c conda-forge -c bioconda homer bedtools samtools bowtie2
conda activate chirp-motif && which findMotifsGenome.pl
```

Step 13 (GO enrichment) needs two more R packages in `chirp-fig`:

```bash
conda install -y -n chirp-fig -c conda-forge -c bioconda bioconductor-clusterprofiler bioconductor-enrichplot
```

---

## 3. Get the scripts onto TSCC

Copy this `tscc` folder into your home directory on TSCC (for example with
`scp`, `git clone`, or by pasting the files). Then go into it:

```bash
cd ~/chirp-seq/src/dubr/tscc   # adjust the path to wherever you put it
```

Step 09 also needs the paper's supplementary tables (DEG lists, differential
H3K27ac peaks). Copy the Excel file from this repo to scratch **once**:

```bash
mkdir -p /tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr/paper
# from your own computer:
scp doc/dubr_paper/Table_S1-S6.xlsx \
    <you>@login.tscc.sdsc.edu:/tscc/lustre/ddn/scratch/<you>/chirp_seq_analysis/dubr/paper/
```

---

## 4. Run the pipeline, one script at a time

Submit each script with `sbatch`. Wait for one to finish before starting the
next (later steps use the output of earlier ones).

```bash
sbatch 00_download.sh          # download the data + genome
sbatch 01_fastqc.sh            # quality check the raw reads
sbatch 01b_trim_fastp.sh       # clean the reads (trim adapters / poly-G)
sbatch 02_align_bowtie2.sh     # align reads to the genome
sbatch 03_filter_merge.sh      # filter reads + combine replicates
sbatch 03b_ucsc_chrom_names.sh # rename "1" -> "chr1" in the merged BAMs (needed for the blacklist)
sbatch 04_peak_calling.sh      # call peaks with MACS2
sbatch 05_confident_peaks.sh   # keep peaks found in BOTH EVEN and ODD
sbatch 06_signal_deeptools.sh  # signal tracks + final peak set (EVEN/ODD ratio) + QC heatmap
```

Steps 00–06 are the ChIRP-seq analysis itself. Steps 07–11 reproduce the
paper's Figures 5, 6 and S8, which overlay the DUBR peaks on the authors'
ChIP-seq / RNA-seq and public ENCODE data (downloaded as ready-made tracks in
step 07 — no re-processing needed). Same rule: one at a time, in order.

```bash
sbatch 07_download_tracks.sh        # GEO + ENCODE bigWigs (~10 GB), replicate averages, refGene
sbatch 08_heatmaps_fig5GH.sh        # Fig 5G (k-means heatmap) + Fig 5H (profiles)
sbatch 09_great_fig5CDE.sh          # Fig 5C, 5D, 5E, EVEN/ODD scatter, S8C   [chirp-fig env]
sbatch 10_tracks_pygenometracks.sh  # Fig 5I, 6A snapshots + DUBR-locus check [chirp-fig env]
                                    #   (tracks are described in tracks_ini/*.ini - edit those to change a figure)
sbatch 11_hes1_fig6DE_S8.sh         # Fig 6D, 6E, S8D, S8J
sbatch 12_motifs_chirp.sh           # probe off-target check + HOMER motifs (~2 h) [chirp-motif env]
sbatch 13_go_enrichment.sh          # GO enrichment (Fig 5F, 6F, S8A, local Metascape) [chirp-fig env]
```

Steps 07 and 09 need internet access from the compute node (GEO/ENCODE
downloads; rGREAT fetches a small TSS table on first use). If a job fails with
a network error, run that step directly on the login node instead
(`bash 07_download_tracks.sh`, or `conda activate chirp-fig && Rscript 09_great_fig5CDE.R`).

### What each step produces

| Script | What it makes |
|---|---|
| `00_download.sh` | Raw FASTQ files + hg19 genome index + blacklist |
| `01_fastqc.sh` | Quality reports (FastQC + one MultiQC summary) |
| `01b_trim_fastp.sh` | Cleaned FASTQ files |
| `02_align_bowtie2.sh` | Aligned BAM files |
| `03_filter_merge.sh` | Filtered + merged BAM files |
| `03b_ucsc_chrom_names.sh` | Same merged BAMs with UCSC chromosome names (`*_merged.chr.bam`) |
| `04_peak_calling.sh` | EVEN and ODD peak files (blacklist-filtered) |
| `05_confident_peaks.sh` | `DUBR_confident_peaks.bed` (EVEN ∩ ODD) + Fig 5C Venn counts |
| `06_signal_deeptools.sh` | bigWig tracks, `DUBR_peaks_final.bed` (EVEN/ODD ratio 0.5–2), QC heatmap |
| `07_download_tracks.sh` | `tracks/avg/*.bw` (histone marks, RNA, authors' ChIRP), `tracks/encode/*` (HES1, GATA-1, …), `ref/annot/` (refGene, TSS) |
| `08_heatmaps_fig5GH.sh` | `figures/fig5G_heatmap_kmeans2.png`, `fig5H_<mark>_profile.png`, `peaks/clusters/C1.bed C2.bed` |
| `09_great_fig5CDE.sh` (+ `.R`) | `fig5C_venn_*.png`, `fig5D_GREAT_*.png/tsv`, `fig5E_venn_*.png`, `fig5E_direct_DUBR_associated_genes.tsv`, `figS8C_*`, Table S4 BEDs |
| `10_tracks_pygenometracks.sh` (+ `tracks_ini/*.ini`) | `figures/tracks/fig5I_<gene>.png`, `fig6A_HES1.png`, `dubr_locus_DUBR_locus.png` |
| `11_hes1_fig6DE_S8.sh` | `fig6D_heatmap_*.png`, `fig6E_donut_*.png` + gene list for Metascape, `figS8D_*`, `figS8J_*` |
| `12_motifs_chirp.sh` (+ `dubr_chirp_probes.fa`) | `motifs/peaks_with_probe_alignment.bed`, `motifs/peaks_with_probe_kmer.txt`, `peaks/DUBR_peaks_final_noProbe.bed`, `motifs/homer_{final,noProbe}_peaks/` (`knownResults.html`, `homerResults.html`), `motifs/homer_C1/` |
| `13_go_enrichment.sh` (+ `.R`) | `go_enrichment/*_GO_BP.tsv`, `*_GO_BP_dotplot.png`, `fig5F_terms_with_TFs.tsv` |

What each figure needs and how to read the results is explained in
`doc/dubr_figures_guide.md` in the main repo.

---

## 5. Check on your jobs

```bash
squeue -u $USER                 # see your jobs (PD = pending, R = running)
scancel <jobid>                 # cancel a job (get the jobid from squeue)
```

When a job finishes, Slurm saves everything it printed to a file named
`<jobname>-<jobid>.out` in this folder (for example `dubr_align-1234567.out`).
Open it to read the log or find errors.

---

## 6. Where the files go (important)

All the big data files are written to the fast scratch disk, **not** your home
folder:

```
/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/dubr/
```

Your home folder is small (100 GB) and is only for code. **Warning:** scratch
files are automatically deleted after 90 days, so once the analysis is done,
copy the important results (the `peaks/` and `figures/` folders) somewhere
permanent.
