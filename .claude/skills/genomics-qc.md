# Quality Control and Preprocessing Skill

You are an expert in quality control and preprocessing of next-generation sequencing data. You help users assess data quality, troubleshoot issues, and prepare data for downstream analysis.

## Quality Control Workflow

```
Raw FASTQ → FastQC → MultiQC → Trimming (if needed) → Post-trim QC → Ready for Analysis
```

## Raw Read Quality Control

### FastQC

```bash
# Single file
fastqc -t 4 -o qc_output/ sample.fastq.gz

# Multiple files
fastqc -t 8 -o qc_output/ *.fastq.gz

# No extraction (just HTML)
fastqc --noextract -t 8 -o qc_output/ *.fastq.gz
```

### FastQC Metrics Interpretation

| Module | Good | Warning | Fail |
|--------|------|---------|------|
| Per base quality | All >Q28 | Some <Q20 | Many <Q20 |
| Per sequence quality | Peak >Q30 | Peak <Q27 | Peak <Q20 |
| Per base N content | <5% | 5-20% | >20% |
| Sequence length | Uniform | Variable | Very variable |
| Overrepresented | <1% | 1-10% | >10% |
| Adapter content | <5% | 5-10% | >10% |
| GC content | Normal curve | Shifted | Bimodal |

### MultiQC (Aggregate Reports)

```bash
# Basic usage
multiqc qc_output/ -o multiqc_report/

# With configuration
multiqc . \
    -o multiqc_report/ \
    --title "Project QC Report" \
    --comment "Sequencing run 2024-01" \
    -f  # Force overwrite

# Specific modules
multiqc . -m fastqc -m star -m featureCounts
```

## Read Trimming

### fastp (Recommended - Fast and Feature-rich)

```bash
# Paired-end
fastp -i R1.fq.gz -I R2.fq.gz \
    -o R1.trimmed.fq.gz -O R2.trimmed.fq.gz \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --thread 8 \
    --html fastp.html \
    --json fastp.json

# Single-end
fastp -i input.fq.gz \
    -o output.fq.gz \
    --detect_adapter_for_se \
    --thread 8 \
    --html fastp.html

# UMI handling
fastp -i R1.fq.gz -I R2.fq.gz \
    -o R1.trimmed.fq.gz -O R2.trimmed.fq.gz \
    --umi \
    --umi_loc read1 \
    --umi_len 8
```

### Trim Galore (Wrapper for Cutadapt)

```bash
# Paired-end
trim_galore --paired --fastqc \
    --cores 4 \
    -o trimmed/ \
    R1.fq.gz R2.fq.gz

# With quality trimming
trim_galore --paired --fastqc \
    --quality 20 \
    --length 36 \
    --cores 4 \
    -o trimmed/ \
    R1.fq.gz R2.fq.gz

# Specify adapter
trim_galore --paired \
    --adapter AGATCGGAAGAGC \
    --adapter2 AGATCGGAAGAGC \
    R1.fq.gz R2.fq.gz
```

### Cutadapt (Direct)

```bash
# Paired-end with Illumina adapters
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20 \
    -m 36 \
    --cores 8 \
    -o R1.trimmed.fq.gz \
    -p R2.trimmed.fq.gz \
    R1.fq.gz R2.fq.gz

# Hard clip from ends
cutadapt -u 10 -u -5 -o output.fq.gz input.fq.gz
```

### BBDuk (BBTools)

```bash
# Quality trimming and adapter removal
bbduk.sh \
    in=R1.fq.gz in2=R2.fq.gz \
    out=R1.trimmed.fq.gz out2=R2.trimmed.fq.gz \
    ref=adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 \
    qtrim=rl trimq=20 \
    minlen=36 \
    threads=8

# Remove PhiX spike-in
bbduk.sh \
    in=input.fq.gz \
    out=clean.fq.gz \
    ref=phix174_ill.ref.fa.gz \
    k=31 hdist=1
```

## Alignment Quality Control

### Samtools Stats

```bash
# Basic stats
samtools flagstat sample.bam > flagstat.txt
samtools stats sample.bam > stats.txt

# Coverage depth
samtools depth -a sample.bam | \
    awk '{sum+=$3} END {print "Average depth:", sum/NR}'

# Quick summary
samtools idxstats sample.bam
```

### Picard Metrics

```bash
# Alignment summary
picard CollectAlignmentSummaryMetrics \
    R=reference.fa \
    I=sample.bam \
    O=alignment_metrics.txt

# Insert size distribution
picard CollectInsertSizeMetrics \
    I=sample.bam \
    O=insert_size_metrics.txt \
    H=insert_size_histogram.pdf

# GC bias
picard CollectGcBiasMetrics \
    I=sample.bam \
    O=gc_bias_metrics.txt \
    CHART=gc_bias_chart.pdf \
    S=gc_summary.txt \
    R=reference.fa

# WGS metrics
picard CollectWgsMetrics \
    I=sample.bam \
    O=wgs_metrics.txt \
    R=reference.fa

# WES/targeted metrics
picard CollectHsMetrics \
    I=sample.bam \
    O=hs_metrics.txt \
    R=reference.fa \
    BAIT_INTERVALS=bait.interval_list \
    TARGET_INTERVALS=target.interval_list

# RNA-seq metrics
picard CollectRnaSeqMetrics \
    I=sample.bam \
    O=rnaseq_metrics.txt \
    REF_FLAT=refFlat.txt \
    STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND
```

### mosdepth (Coverage Analysis)

```bash
# Genome-wide coverage
mosdepth -t 4 -n --fast-mode sample sample.bam

# By chromosome
mosdepth -t 4 -b regions.bed sample sample.bam

# Per-base coverage (for small regions)
mosdepth -t 4 -b target.bed --by 1 sample sample.bam

# Quantize coverage into bins
mosdepth -t 4 -n --quantize 0:1:5:10:50:100:500 sample sample.bam
```

### Qualimap

```bash
# BAM QC
qualimap bamqc \
    -bam sample.bam \
    -outdir qualimap_output/ \
    -nt 8 \
    --java-mem-size=16G

# RNA-seq QC
qualimap rnaseq \
    -bam sample.bam \
    -gtf genes.gtf \
    -outdir qualimap_rnaseq/ \
    --java-mem-size=16G

# Multi-BAM QC
qualimap multi-bamqc \
    -d sample_list.txt \
    -outdir qualimap_multi/
```

## Sample Quality Assessment

### VerifyBamID (Contamination Check)

```bash
# Check contamination
verifybamid2 \
    --Ref reference.fa \
    --BamFile sample.bam \
    --SVDPrefix resource_prefix \
    --Output sample.verifybam

# Interpret: FREEMIX > 0.03 suggests contamination
```

### Somalier (Relatedness and Sex Check)

```bash
# Extract sites
somalier extract -d extracted/ \
    --sites sites.vcf.gz \
    -f reference.fa \
    *.bam

# Relate samples
somalier relate -o results extracted/*.somalier

# Outputs ancestry, relatedness, sex inference
```

### Peddy (Pedigree/Sex Check)

```bash
peddy -p 4 \
    --plot \
    -o peddy_output \
    joint.vcf.gz \
    sample.ped
```

## Variant Quality Control

### VCF Statistics

```bash
# bcftools stats
bcftools stats -s - annotated.vcf.gz > vcf_stats.txt

# Per-sample stats
bcftools stats -s sample1,sample2 joint.vcf.gz > sample_stats.txt

# Ti/Tv ratio
bcftools stats annotated.vcf.gz | grep "TSTV"

# Plot VCF stats
plot-vcfstats -p vcf_plots/ vcf_stats.txt
```

### GATK Metrics

```bash
# Variant evaluation
gatk VariantEval \
    -R reference.fa \
    --eval variants.vcf.gz \
    --dbsnp dbsnp.vcf.gz \
    -O eval_metrics.grp

# Count variants
gatk CountVariants \
    -V variants.vcf.gz
```

## Batch QC and Outlier Detection

### Principal Component Analysis

```r
library(ggplot2)
library(dplyr)

# Load QC metrics
metrics <- read.delim("combined_metrics.txt")

# PCA on QC metrics
qc_pca <- prcomp(metrics %>% select_if(is.numeric), scale = TRUE)

# Plot
ggplot(as.data.frame(qc_pca$x), aes(PC1, PC2, color = metrics$batch)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "QC PCA - Batch Effect Check")

# Identify outliers (>3 SD from mean)
outliers <- metrics %>%
  mutate(across(where(is.numeric), ~ abs(scale(.)) > 3)) %>%
  filter(if_any(everything(), ~ . == TRUE))
```

### Automated QC with MultiQC

```yaml
# multiqc_config.yaml
title: "Project QC Report"
show_analysis_paths: False
show_analysis_time: False

# Highlight thresholds
fastqc_config:
  fastqc_theoretical_gc: "hg38_genome"

# Custom thresholds
table_columns_visible:
  FastQC:
    percent_duplicates: True
    percent_gc: True
    avg_sequence_length: True
    percent_fails: True

custom_data:
  my_genstats:
    plot_type: 'generalstats'
    pconfig:
      - sample_qualdepth:
          title: 'Quality Depth'
          min: 0
          max: 100
```

## Common QC Issues and Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| Low quality scores at 3' end | Normal sequencing degradation | Trim low-quality bases |
| High adapter content | Short inserts / over-sequencing | Adapter trimming |
| High duplication | PCR over-amplification | Mark/remove duplicates |
| GC bias | Library prep issue | Use GC-aware tools |
| Contamination | Sample mix-up / cross-contamination | Verify with contamination tools |
| Wrong sex | Sample swap | Check with sex inference tools |
| Batch effects | Technical variation | Include batch in analysis |
| Low mapping rate | Wrong reference / contamination | Check species, adapter content |
| High rRNA (RNA-seq) | Incomplete rRNA depletion | Filter in silico |

## QC Report Template

```markdown
# Sequencing QC Report

## Summary
- Total samples: X
- Passed QC: Y
- Failed QC: Z

## Key Metrics

### Read Quality
| Sample | Total Reads | Q30 (%) | Adapter (%) | Status |
|--------|-------------|---------|-------------|--------|
| S1     | 50M         | 95      | 2           | PASS   |
| S2     | 45M         | 88      | 15          | WARN   |

### Alignment
| Sample | Mapping (%) | Duplicates (%) | Coverage |
|--------|-------------|----------------|----------|
| S1     | 98          | 15             | 35x      |
| S2     | 92          | 45             | 28x      |

## Failed Samples
- S3: High contamination (FREEMIX=0.08)
- S4: Sex mismatch (expected F, observed M)

## Recommendations
1. Re-sequence S2 due to high adapter content
2. Investigate S3 and S4 sample swaps
```

## Automated QC Pipeline (Nextflow)

```bash
# nf-core/fastqc
nextflow run nf-core/fastqc -profile docker --input samplesheet.csv

# For RNA-seq QC
nextflow run nf-core/rnaseq -profile docker \
    --input samplesheet.csv \
    --genome GRCh38 \
    --skip_alignment  # QC only
```
