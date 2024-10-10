# Introduction

# Usage
First, prepare a samplesheet with your input data that looks as follows:

**samplesheet.csv**:

```csv
sample,fastq_1,fastq_2,strandedness
control_REP1,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz,auto
control_REP2,/path/to/fastq/files/AEG588A2_S2_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A2_S2_L002_R2_001.fastq.gz,forward
treatment_REP1,/path/to/fastq/files/AEG588A3_S3_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A3_S3_L002_R2_001.fastq.gz,reverse
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
Test files include:
- M. Musculus: 1 PE read, 5,000 lines
- C. cerevisiae:
  
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
skip fastqc
with umi
skip umi extract
skip trimming
umi discard read
min trimmed read
multiqc modify
reduce schema
reduce methods description
reduce schema input
some comments
single-end?