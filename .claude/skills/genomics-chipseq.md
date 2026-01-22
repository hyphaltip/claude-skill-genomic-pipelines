# ChIP-seq and ATAC-seq Analysis Pipeline Skill

You are an expert in chromatin immunoprecipitation sequencing (ChIP-seq) and Assay for Transposase-Accessible Chromatin sequencing (ATAC-seq) analysis. You help users design and implement epigenomic analysis workflows.

## Workflow Overview

### ChIP-seq
```
FASTQ → QC → Alignment → Filtering → Peak Calling → Annotation → Differential Binding → Motif Analysis
```

### ATAC-seq
```
FASTQ → QC → Alignment → Filtering → Peak Calling → Nucleosome Positioning → Footprinting → TF Analysis
```

## ChIP-seq Pipeline

### 1. Quality Control and Preprocessing

```bash
# FastQC
fastqc -t 4 -o qc/ *.fastq.gz

# Adapter trimming
trim_galore --fastqc -o trimmed/ sample.fastq.gz

# For paired-end
trim_galore --paired --fastqc -o trimmed/ sample_R1.fastq.gz sample_R2.fastq.gz
```

### 2. Alignment

```bash
# Build Bowtie2 index
bowtie2-build genome.fa genome_index

# Align single-end
bowtie2 -p 8 \
    -x genome_index \
    -U sample.fastq.gz \
    --very-sensitive \
    2> sample.bowtie2.log | \
    samtools view -bS -q 30 - | \
    samtools sort -@ 4 -o sample.sorted.bam -

# Align paired-end
bowtie2 -p 8 \
    -x genome_index \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    -X 1000 \
    2> sample.bowtie2.log | \
    samtools view -bS -F 4 -q 30 - | \
    samtools sort -@ 4 -o sample.sorted.bam -

# Index
samtools index sample.sorted.bam
```

### 3. Post-Alignment Filtering

```bash
# Remove duplicates
picard MarkDuplicates \
    INPUT=sample.sorted.bam \
    OUTPUT=sample.dedup.bam \
    METRICS_FILE=sample.dup_metrics.txt \
    REMOVE_DUPLICATES=true

# Or use samtools
samtools markdup -r sample.sorted.bam sample.dedup.bam

# Remove blacklist regions
bedtools intersect -v -abam sample.dedup.bam -b blacklist.bed > sample.filtered.bam

# For paired-end: keep properly paired reads
samtools view -f 2 -b sample.dedup.bam > sample.proper.bam
```

### 4. Quality Metrics

```bash
# Alignment stats
samtools flagstat sample.filtered.bam

# Fragment size distribution (paired-end)
picard CollectInsertSizeMetrics \
    INPUT=sample.filtered.bam \
    OUTPUT=insert_size_metrics.txt \
    HISTOGRAM_FILE=insert_size_histogram.pdf

# Fingerprint plot (for input comparison)
plotFingerprint -b sample.filtered.bam input.filtered.bam \
    --labels Sample Input \
    -o fingerprint.png

# Cross-correlation (phantompeakqualtools)
run_spp.R -c=sample.filtered.bam -savp -out=spp_metrics.txt
```

### 5. Peak Calling with MACS2

```bash
# Narrow peaks (TF ChIP-seq)
macs2 callpeak \
    -t sample.filtered.bam \
    -c input.filtered.bam \
    -f BAM \
    -g hs \
    -n sample \
    -B \
    --SPMR \
    -q 0.05 \
    --outdir peaks/

# Broad peaks (histone marks: H3K27me3, H3K36me3, H3K9me3)
macs2 callpeak \
    -t sample.filtered.bam \
    -c input.filtered.bam \
    -f BAM \
    -g hs \
    -n sample \
    --broad \
    --broad-cutoff 0.1 \
    -B \
    --SPMR \
    --outdir peaks/

# Paired-end
macs2 callpeak \
    -t sample.filtered.bam \
    -c input.filtered.bam \
    -f BAMPE \
    -g hs \
    -n sample \
    --outdir peaks/
```

### 6. Peak Annotation (R)

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read peaks
peaks <- readPeakFile("sample_peaks.narrowPeak")

# Annotate peaks
peakAnno <- annotatePeak(peaks, tssRegion = c(-3000, 3000),
                         TxDb = txdb, annoDb = "org.Hs.eg.db")

# Visualize
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)

# GO enrichment of target genes
genes <- as.data.frame(peakAnno)$geneId
ego <- enrichGO(gene = genes,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

# Profile around TSS
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peaks, windows = promoter)
plotAvgProf(tagMatrix, xlim = c(-3000, 3000))
```

### 7. Differential Binding Analysis (R)

```r
library(DiffBind)

# Create sample sheet
samples <- read.csv("samplesheet.csv")
# Columns: SampleID, Tissue, Factor, Condition, Replicate, bamReads, ControlID, bamControl, Peaks

# Create DBA object
dba <- dba(sampleSheet = samples)

# Count reads in peaks
dba <- dba.count(dba)

# Normalize
dba <- dba.normalize(dba)

# Contrast
dba <- dba.contrast(dba, categories = DBA_CONDITION)

# Differential analysis
dba <- dba.analyze(dba)

# Get results
results <- dba.report(dba)

# Visualization
dba.plotMA(dba)
dba.plotVolcano(dba)
dba.plotHeatmap(dba)
```

### 8. Motif Analysis

```bash
# HOMER findMotifsGenome
findMotifsGenome.pl peaks.bed hg38 motif_output/ -size 200 -mask

# MEME-ChIP
meme-chip -oc meme_output \
    -db motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant.meme \
    peak_sequences.fa
```

## ATAC-seq Pipeline

### 1. Alignment (with proper flags for ATAC)

```bash
# Align with Bowtie2
bowtie2 -p 8 \
    -x genome_index \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    -X 2000 \
    -I 10 \
    2> sample.bowtie2.log | \
    samtools view -bS -F 4 -q 30 - | \
    samtools sort -@ 4 -o sample.sorted.bam -
```

### 2. ATAC-specific Filtering

```bash
# Remove duplicates
picard MarkDuplicates \
    INPUT=sample.sorted.bam \
    OUTPUT=sample.dedup.bam \
    METRICS_FILE=dup_metrics.txt \
    REMOVE_DUPLICATES=true

# Remove mitochondrial reads
samtools view -h sample.dedup.bam | \
    grep -v chrM | \
    samtools view -b - > sample.noMT.bam

# Remove blacklist regions
bedtools intersect -v -abam sample.noMT.bam -b blacklist.bed > sample.filtered.bam

# Shift reads for Tn5 offset (+4 on + strand, -5 on - strand)
alignmentSieve -b sample.filtered.bam \
    -o sample.shifted.bam \
    --ATACshift

# Or with samtools (manual shift)
samtools view -h sample.filtered.bam | \
    awk 'BEGIN {OFS="\t"} {
        if ($1 ~ /^@/) {print; next}
        if (and($2, 16)) {$4 = $4 - 5}
        else {$4 = $4 + 4}
        print
    }' | samtools view -b - > sample.shifted.bam
```

### 3. Peak Calling for ATAC-seq

```bash
# MACS2 for ATAC-seq (no control)
macs2 callpeak \
    -t sample.shifted.bam \
    -f BAMPE \
    -g hs \
    -n sample \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    -B \
    --SPMR \
    --keep-dup all \
    -q 0.05 \
    --outdir peaks/

# HMMRATAC (specific for ATAC-seq, better at nucleosome detection)
HMMRATAC -b sample.filtered.bam \
    -i sample.filtered.bam.bai \
    -g genome.info \
    -o sample_hmmratac
```

### 4. Nucleosome-Free Region Analysis

```bash
# Separate reads by fragment size
samtools view -h sample.filtered.bam | \
    awk 'substr($0,1,1)=="@" || ($9>=0 && $9<=100) || ($9<=0 && $9>=-100)' | \
    samtools view -b - > nucleosome_free.bam

samtools view -h sample.filtered.bam | \
    awk 'substr($0,1,1)=="@" || ($9>=180 && $9<=247) || ($9<=-180 && $9>=-247)' | \
    samtools view -b - > mononucleosome.bam
```

### 5. TSS Enrichment Score

```r
library(ATACseqQC)

# Calculate TSS enrichment
tsse <- TSSEscore(sample.bam, txdb)
print(paste("TSS Enrichment Score:", tsse$TSSEscore))

# Fragment size distribution
fragSizeDist(sample.bam, "Sample")
```

### 6. Footprinting Analysis

```bash
# TOBIAS for TF footprinting
# Correct Tn5 bias
TOBIAS ATACorrect \
    --bam sample.filtered.bam \
    --genome genome.fa \
    --peaks peaks.bed \
    --outdir tobias/ \
    --cores 8

# Score footprints
TOBIAS ScoreBigwig \
    --signal tobias/sample_corrected.bw \
    --regions peaks.bed \
    --output tobias/sample_footprints.bw \
    --cores 8

# Differential footprinting
TOBIAS BINDetect \
    --motifs JASPAR2022_CORE_vertebrates.meme \
    --signals condition1_footprints.bw condition2_footprints.bw \
    --genome genome.fa \
    --peaks peaks.bed \
    --outdir tobias/bindetect \
    --cores 8
```

## nf-core Pipelines

### ChIP-seq
```bash
nextflow run nf-core/chipseq \
    -profile docker \
    --input samplesheet.csv \
    --genome GRCh38 \
    --narrow_peak  # or --broad_peak
```

### ATAC-seq
```bash
nextflow run nf-core/atacseq \
    -profile docker \
    --input samplesheet.csv \
    --genome GRCh38
```

## Quality Metrics

### ChIP-seq
| Metric | Good | Concern |
|--------|------|---------|
| Total reads | >20M | <10M |
| Mapping rate | >80% | <70% |
| Duplicate rate | <20% | >40% |
| NSC (normalized strand coefficient) | >1.1 | <1.05 |
| RSC (relative strand correlation) | >1.0 | <0.8 |
| FRiP (fraction reads in peaks) | >1% | <1% |

### ATAC-seq
| Metric | Good | Concern |
|--------|------|---------|
| Total reads | >50M | <25M |
| Mapping rate | >80% | <70% |
| Mitochondrial rate | <20% | >50% |
| Duplicate rate | <30% | >50% |
| TSS enrichment | >6 | <4 |
| FRiP | >20% | <10% |
| NFR/mono ratio | >2.5 | <1.5 |

## Common Issues

1. **Low FRiP**: Weak immunoprecipitation, antibody issues
2. **No peaks**: Check input control quality, signal-to-noise
3. **High background**: Increase sequencing depth, better antibody
4. **Batch effects**: Use spike-in controls (e.g., Drosophila chromatin)
5. **High MT reads (ATAC)**: Cell death during prep, improve protocol
