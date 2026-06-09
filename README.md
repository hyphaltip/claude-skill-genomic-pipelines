# Genomics Skills for Claude Code

A collection of bioinformatics skills for [Claude Code](https://claude.ai/claude-code) that provide expert guidance on genomic data analysis pipelines.

## Installation

Clone this repository and Claude Code will automatically detect the skills when running in this directory:

```bash
git clone https://github.com/hyphaltip/claude_skill_genomics.git
cd claude_skill_genomics
```

Or copy the `.claude` folder to your existing project:

```bash
cp -r claude_skill_genomics/.claude /path/to/your/project/
```

## Available Skills

| Skill | Command | Description |
|-------|---------|-------------|
| Genomics | `/genomics` | Main skill - overview, tool selection, routing |
| WGS/WES | `/genomics:wgs` | Variant calling (GATK, DeepVariant, Mutect2) |
| RNA-seq | `/genomics:rnaseq` | Bulk & single-cell RNA-seq (DESeq2, Seurat, Scanpy) |
| ChIP/ATAC-seq | `/genomics:chipseq` | Peak calling, footprinting, motif analysis |
| Annotation | `/genomics:annotation` | Variant annotation (VEP, SnpEff, ANNOVAR) |
| QC | `/genomics:qc` | Quality control & preprocessing (FastQC, MultiQC) |
| CNV/SV | `/genomics:cnv` | Copy number & structural variants (GATK CNV, CNVkit, Manta) |

## Usage

In any Claude Code session:

```
/genomics           # General bioinformatics help
/genomics:wgs       # WGS/WES variant calling guidance
/genomics:rnaseq    # RNA-seq analysis help
/genomics:cnv       # CNV/SV detection pipelines
```

## What's Included

Each skill provides:

- Complete command-line examples with realistic parameters
- Tool comparison tables for selecting the right tool
- nf-core pipeline recommendations for production use
- Quality thresholds and troubleshooting guides
- R/Python code for downstream analysis
- Best practices following GATK, ENCODE, and community standards

## Covered Tools

### Alignment
BWA-MEM2, Bowtie2, STAR, HISAT2, Minimap2

### Variant Calling
GATK HaplotypeCaller, DeepVariant, Mutect2, Strelka2, FreeBayes, bcftools

### CNV/SV
GATK CNV, CNVkit, Manta, DELLY, GRIDSS, LUMPY

### RNA-seq
Salmon, kallisto, STAR, featureCounts, DESeq2, edgeR, Seurat, Scanpy

### ChIP/ATAC-seq
MACS2, HOMER, DiffBind, TOBIAS

### Annotation
VEP, SnpEff, ANNOVAR, AnnotSV

### QC
FastQC, MultiQC, fastp, Picard, mosdepth, Qualimap

## File Structure

```
.claude/
├── settings.json          # Skill registration
└── skills/
    ├── genomics.md        # Main skill
    ├── genomics-wgs.md    # WGS/WES pipelines
    ├── genomics-rnaseq.md # RNA-seq analysis
    ├── genomics-chipseq.md# ChIP/ATAC-seq
    ├── genomics-annotation.md # Variant annotation
    ├── genomics-qc.md     # Quality control
    └── genomics-cnv.md    # CNV/SV analysis
```

## Requirements

- [Claude Code CLI](https://claude.ai/claude-code)
- No additional dependencies (skills are prompt-based)

## Contributing

Contributions welcome! To add or improve skills:

1. Fork the repository
2. Edit or add skill files in `.claude/skills/`
3. Register new skills in `.claude/settings.json`
4. Submit a pull request

## License

MIT
