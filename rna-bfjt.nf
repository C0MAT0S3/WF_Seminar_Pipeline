// Import modules
include { STAR_ALIGN } from './modules/nf-core/star/align/main' 
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main'

// Import subworkflows
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from './subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { ALIGN_STAR } from './subworkflows/local/align_star/main'
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome/main'
include { fromSamplesheet                  } from 'plugin/nf-validation'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = getGenomeAttribute('fasta')
params.additional_fasta = getGenomeAttribute('additional_fasta')
params.transcript_fasta = getGenomeAttribute('transcript_fasta')
// params.gff              = getGenomeAttribute('gff')
params.gtf              = getGenomeAttribute('gtf')
// params.gene_bed         = getGenomeAttribute('bed12')
// params.bbsplit_index    = getGenomeAttribute('bbsplit')
// params.sortmerna_index  = getGenomeAttribute('sortmerna')
params.star_index       = getGenomeAttribute('star')
// params.rsem_index       = getGenomeAttribute('rsem')
// params.hisat2_index     = getGenomeAttribute('hisat2')
// params.salmon_index     = getGenomeAttribute('salmon')
// params.kallisto_index   = getGenomeAttribute('kallisto')

workflow {
// Reading the samplesheet file
Channel
    .fromSamplesheet("input")
    .map {
        meta, fastq_1, fastq_2 ->
            if (!fastq_2) {
                return [ meta + [ single_end:true ], [ fastq_1 ] ]
            } else {
                return [ meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
            }
    }
    .set{ reads_ch }

// QC
FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
    reads_ch, // reads             // channel: [ val(meta), [ reads ] ]
    false, // skip_fastqc       // boolean: true/false
    false, // with_umi          // boolean: true/false
    true, // skip_umi_extract  // boolean: true/false
    false, // skip_trimming     // boolean: true/false
    0, // umi_discard_read  // integer: 0, 1 or 2
    1 // min_trimmed_reads // integer: > 0
)

// Prepare Genome
// TODO fasta empty
PREPARE_GENOME(
    params.fasta,
    params.gtf,
    params.star_index,
    params.gencode
)
PREPARE_GENOME.out.star_index.view()

// Aligning
ALIGN_STAR(
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads, // reads               // channel: [ val(meta), [ reads ] ]
    PREPARE_GENOME.out.star_index.map { [ [:], it ] },               // channel: [ val(meta), [ index ] ]
    PREPARE_GENOME.out.gtf.map { [ [:], it ] },                // channel: [ val(meta), [ gtf ] ]
    false, // star_ignore_sjdbgtf // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    '', // seq_platform        // string : sequencing platform
    '', // seq_center          // string : sequencing center
    true, // is_aws_igenome      // boolean: whether the genome files are from AWS iGenomes
    PREPARE_GENOME.out.fasta               // channel: /path/to/fasta
)

// MultiQC report (optional)
//TODO
// Markduplicates (optional)
//TODO
}

//
// Get attribute from genome config file e.g. fasta
//
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}