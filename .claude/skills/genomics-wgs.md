# WGS/WES Variant Calling Pipeline Skill

You are an expert in whole genome sequencing (WGS) and whole exome sequencing (WES) analysis pipelines. You help users design and implement variant calling workflows following GATK best practices.

## Workflow Overview

```
FASTQ → QC → Alignment → Post-processing → Variant Calling → Filtering → Annotation
```

## Standard WGS/WES Pipeline

### 1. Quality Control (Raw Reads)

```bash
# FastQC on raw reads
fastqc -t 4 -o qc_raw/ sample_R1.fastq.gz sample_R2.fastq.gz

# Optional: Adapter trimming (usually not needed for modern data)
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
    -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
    --detect_adapter_for_pe \
    --thread 4 \
    --html fastp_report.html
```

### 2. Alignment

```bash
# BWA-MEM2 alignment (faster than BWA-MEM)
bwa-mem2 mem -t 16 \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
    ${REFERENCE} \
    ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz | \
    samtools sort -@ 8 -m 2G -o ${SAMPLE}.sorted.bam -

# Index BAM
samtools index -@ 4 ${SAMPLE}.sorted.bam
```

### 3. Post-Alignment Processing

```bash
# Mark duplicates
gatk MarkDuplicates \
    -I ${SAMPLE}.sorted.bam \
    -O ${SAMPLE}.dedup.bam \
    -M ${SAMPLE}.dup_metrics.txt \
    --CREATE_INDEX true

# Base Quality Score Recalibration (BQSR)
# Step 1: Generate recalibration table
gatk BaseRecalibrator \
    -R ${REFERENCE} \
    -I ${SAMPLE}.dedup.bam \
    --known-sites ${DBSNP} \
    --known-sites ${MILLS_INDELS} \
    --known-sites ${KNOWN_INDELS} \
    -O ${SAMPLE}.recal_data.table

# Step 2: Apply recalibration
gatk ApplyBQSR \
    -R ${REFERENCE} \
    -I ${SAMPLE}.dedup.bam \
    --bqsr-recal-file ${SAMPLE}.recal_data.table \
    -O ${SAMPLE}.bqsr.bam
```

### 4. Variant Calling

#### Option A: GATK HaplotypeCaller (Standard)

```bash
# Single sample calling
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I ${SAMPLE}.bqsr.bam \
    -O ${SAMPLE}.g.vcf.gz \
    -ERC GVCF \
    --native-pair-hmm-threads 4

# Joint genotyping (for cohorts)
gatk GenomicsDBImport \
    -V sample1.g.vcf.gz \
    -V sample2.g.vcf.gz \
    --genomicsdb-workspace-path genomicsdb \
    -L intervals.bed

gatk GenotypeGVCFs \
    -R ${REFERENCE} \
    -V gendb://genomicsdb \
    -O cohort.vcf.gz
```

#### Option B: DeepVariant (ML-based, highly accurate)

```bash
# Using Docker
docker run -v "${PWD}":"/input" -v "${PWD}":"/output" \
    google/deepvariant:1.6.0 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/reference.fa \
    --reads=/input/${SAMPLE}.bqsr.bam \
    --output_vcf=/output/${SAMPLE}.vcf.gz \
    --output_gvcf=/output/${SAMPLE}.g.vcf.gz \
    --num_shards=16
```

### 5. Variant Filtering

#### GATK VQSR (for large cohorts, 30+ samples)

```bash
# SNPs
gatk VariantRecalibrator \
    -R ${REFERENCE} \
    -V cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O snp.recal \
    --tranches-file snp.tranches

gatk ApplyVQSR \
    -R ${REFERENCE} \
    -V cohort.vcf.gz \
    --recal-file snp.recal \
    --tranches-file snp.tranches \
    -mode SNP \
    --truth-sensitivity-filter-level 99.5 \
    -O cohort.snp_recal.vcf.gz

# INDELs
gatk VariantRecalibrator \
    -R ${REFERENCE} \
    -V cohort.snp_recal.vcf.gz \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 mills.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    -O indel.recal \
    --tranches-file indel.tranches

gatk ApplyVQSR \
    -R ${REFERENCE} \
    -V cohort.snp_recal.vcf.gz \
    --recal-file indel.recal \
    --tranches-file indel.tranches \
    -mode INDEL \
    --truth-sensitivity-filter-level 99.0 \
    -O cohort.filtered.vcf.gz
```

#### Hard Filtering (for small sample sets)

```bash
# SNPs
gatk SelectVariants -V raw.vcf.gz -select-type SNP -O snps.vcf.gz

gatk VariantFiltration \
    -V snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O snps.filtered.vcf.gz

# INDELs
gatk SelectVariants -V raw.vcf.gz -select-type INDEL -O indels.vcf.gz

gatk VariantFiltration \
    -V indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O indels.filtered.vcf.gz

# Merge
gatk MergeVcfs \
    -I snps.filtered.vcf.gz \
    -I indels.filtered.vcf.gz \
    -O filtered.vcf.gz
```

## Reference Files (GRCh38)

```bash
# Option 1: Broad Institute GATK Resource Bundle (public FTP - no auth required)
# https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references
BROAD_FTP="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"

wget ${BROAD_FTP}/Homo_sapiens_assembly38.fasta
wget ${BROAD_FTP}/Homo_sapiens_assembly38.fasta.fai
wget ${BROAD_FTP}/Homo_sapiens_assembly38.dict
wget ${BROAD_FTP}/Homo_sapiens_assembly38.dbsnp138.vcf
wget ${BROAD_FTP}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget ${BROAD_FTP}/Homo_sapiens_assembly38.known_indels.vcf.gz

# Option 2: NCBI (always public)
# Reference genome (GRCh38 analysis set)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# dbSNP (latest)
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

# Option 3: Ensembl (public)
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/variation/vcf/homo_sapiens/homo_sapiens-chr*.vcf.gz

# Option 4: AWS iGenomes (public S3, no auth)
aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ . --recursive

# Option 5: Use nf-core (auto-downloads all references)
nextflow run nf-core/sarek --genome GATK.GRCh38 ...
```

## Quality Metrics to Check

| Metric | Good Value | Concern |
|--------|------------|---------|
| Mapping rate | >95% | <90% |
| Duplicate rate | <20% (WGS), <30% (WES) | >40% |
| Mean coverage | 30x (WGS), 100x (WES) | <20x |
| Coverage uniformity | >80% at 20x | <70% |
| Insert size | 300-500bp | Bimodal |
| Ti/Tv ratio | ~2.1 (WGS), ~3.0 (WES) | <1.8 |

## nf-core/sarek Pipeline

For production use, recommend nf-core/sarek:

```bash
nextflow run nf-core/sarek \
    -profile docker \
    --input samplesheet.csv \
    --genome GRCh38 \
    --tools haplotypecaller,deepvariant,snpeff,vep \
    --wes  # Add for WES data
```

Samplesheet format:
```csv
patient,sample,lane,fastq_1,fastq_2
patient1,sample1,lane1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
```

## Somatic Variant Calling (Tumor/Normal)

```bash
# Mutect2 for somatic variants
gatk Mutect2 \
    -R ${REFERENCE} \
    -I tumor.bam \
    -I normal.bam \
    -normal normal_sample_name \
    --germline-resource af-only-gnomad.vcf.gz \
    --panel-of-normals pon.vcf.gz \
    -O somatic.vcf.gz

# Filter somatic calls
gatk FilterMutectCalls \
    -R ${REFERENCE} \
    -V somatic.vcf.gz \
    -O somatic.filtered.vcf.gz
```

## Common Issues & Solutions

1. **Low mapping rate**: Check read quality, contamination, wrong reference
2. **High duplicate rate**: PCR over-amplification, low library complexity
3. **Uneven coverage**: GC bias, capture issues (WES)
4. **Wrong Ti/Tv**: Possible contamination or systematic errors
5. **Missing variants**: Check coverage at locus, filtering too strict
