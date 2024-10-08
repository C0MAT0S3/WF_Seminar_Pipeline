# Introduction

# Usage
First, prepare a samplesheet with your input data that looks as follows:

**samplesheet.csv**:

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,L002_R1_001.fastq.gz,L002_R2_001.fastq.gz,auto
CONTROL_REP1,L003_R1_001.fastq.gz,L003_R2_001.fastq.gz,forward
CONTROL_REP1,L004_R1_001.fastq.gz,L004_R2_001.fastq.gz,reverse
```

Now, you can run the pipeline using:

```bash
nextflow run rna-bfjt \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --genome GRCh37
```

# Pipeline Output

# Example

# Contributions

# TODO
1. Specify samplesheet for our pipeline
2. Quality Control before trimming (FASTQC)
3. Trimming (TRIMGALORE)
4. Quality Control after trimming (FASTQC)
5. Download Ref. genome (NCBIGENOMEDOWNLOAD)
6. Aligning to Ref. genome (STAR)
7. (MULTIQC for aggregated results)
8. (Mark duplicates)