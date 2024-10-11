# Introduction

# Usage
1. First, prepare a samplesheet with your input read data that looks as follows:

**samplesheet.csv**:

```csv
sample,fastq_1,fastq_2,strandedness
control_REP1,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz,auto
control_REP2,/path/to/fastq/files/AEG588A2_S2_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A2_S2_L002_R2_001.fastq.gz,forward
treatment_REP1,/path/to/fastq/files/AEG588A3_S3_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A3_S3_L002_R2_001.fastq.gz,reverse
```
1. Make sure either Docker or Singularity is set up to be able to pull containers. You can verify their installation by running one of these commands:
```bash
docker run hello-world
singularity run library://sylabsed/examples/lolcow
```
1. Now, you can run the pipeline using:
```bash
nextflow run ./rna-bfjt.nf \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --genome <GENOME> \
    --profile <docker/singularity>
```

## Required Parameters
- **input**: Path to the input samplesheet file (e.g., `.csv` format).
- **outdir**: Output directory where results will be saved.
- **genome**: Reference genome for alignment.
- **profile**: Nextflow profile to define execution environment (e.g., `docker`, `singularity`).

## Optional Parameters
### Maximum Computing Resources
- **max_memory**: Maximum memory allocated for each process. Default: `12.GB`
- **max_cpus**: Maximum CPU cores allocated for each process. Default: `6`
- **max_time**: Maximum execution time for each process. Default: `12.h`

### Boilerplate Options
- **publish_dir_mode**: Method for publishing files to the output directory. Default: `copy`

### References
- **splicesites**: Path to the splice sites file required for some aligners.
- **save_reference**: Option to save the reference genome index after creation. Default: `false`
- **igenomes_base**: Base URL or path to iGenomes reference files. Default: `s3://ngi-igenomes/igenomes/`

### Quality Control
- **skip_fastqc**: Option to skip FastQC analysis. Default: `false`

### UMI Handling
- **with_umi**: Enable UMI-based deduplication. Default: `false`
- **skip_umi_extract**: Skip UMI extraction step. Default: `false`
- **umitools_extract_method**: UMI extraction method, either `string` or `regex`. Default: `string`
- **umitools_grouping_method**: Method for UMI grouping, options include `unique`, `percentile`, `cluster`, `adjacency`, or `directional`. Default: `directional`
- **umitools_dedup_stats**: Generate deduplication statistics. Default: `false`
- **umitools_bc_pattern**: UMI barcode pattern.
- **umitools_bc_pattern2**: UMI barcode pattern for read 2.
- **umitools_umi_separator**: Separator character for UMI in read names.
- **umi_discard_read**: Discard read 1 or read 2 based on UMI extraction.
- **save_umi_intermeds**: Save UMI intermediate files. Default: `false`

### Trimming
- **trimmer**: Tool for trimming, e.g., `trimgalore`. Default: `trimgalore`
- **min_trimmed_reads**: Minimum number of reads to retain after trimming. Default: `1000`
- **extra_trimgalore_args**: Additional arguments for the trimming tool.
- **save_trimmed**: Save trimmed reads. Default: `false`
- **skip_trimming**: Option to skip the trimming step. Default: `false`

### Alignment
- **aligner**: Alignment tool, e.g. `hisat2`. Default: `hisat2`
- **seq_center**: Sequencing center information for BAM files.
- **hisat2_build_memory**: Memory allocation for HISAT2 index build. Default: `12GB`
- **save_unaligned**: Save unaligned reads. Default: `false`
- **save_align_intermeds**: Save intermediate alignment files. Default: `false`
- **skip_markduplicates**: Skip the MarkDuplicates step. Default: `false`

### MultiQC Options
- **multiqc_config**: Path to a custom MultiQC configuration file.
- **multiqc_logo**: Path to a logo file for MultiQC.
- **multiqc_methods_description**: Path to a custom MultiQC methods description file.

# Pipeline Output
The results of the analyses are aggregated using MultiQC and can be found under `<OUTDIR>/multiqc`. Intermediate results can be found under `<OUTDIR>` in their respective folders.

# Example
2 test sets are included in `examples` with a template samplesheet:

| Organism       | SRA ID        | Read Count |
|----------------|---------------|------------|
| M. musculus    | SRR23195511   |  1,250     |
| S. cerevisiae  | SRR6357070    | 50,000     |
| S. cerevisiae  | SRR6357071    | 50,000     |
| S. cerevisiae  | SRR6357072    | 50,000     |
| S. cerevisiae  | SRR6357073    | 50,000     |
| S. cerevisiae  | SRR6357074    | 50,000     |
| S. cerevisiae  | SRR6357075    | 50,000     |
| S. cerevisiae  | SRR6357076    | 50,000     |
| S. cerevisiae  | SRR6357077    | 50,000     |
| S. cerevisiae  | SRR6357078    | 50,000     |
| S. cerevisiae  | SRR6357079    | 50,000     |
| S. cerevisiae  | SRR6357080    | 50,000     |
| S. cerevisiae  | SRR6357081    | 50,000     |

# Contributions
The pipeline was built by:
- Bhavya Kandi, 
- Franziska Sassen, 
- Jan Steiger, 
- Thomas Vogel