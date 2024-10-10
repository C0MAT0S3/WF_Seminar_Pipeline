// Import methods
include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet                  } from 'plugin/nf-validation'
include { multiqcTsvFromList             } from './subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { paramsSummaryMultiqc             } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText         } from './subworkflows/local/utils_nfcore_methods'

// Import modules
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'    
include { MULTIQC } from './modules/local/multiqc/main'                                 

// Import subworkflows
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from './subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome/main'
include { FASTQ_ALIGN_HISAT2 } from './subworkflows/nf-core/fastq_align_hisat2/main'
include { BAM_MARKDUPLICATES_PICARD        } from './subworkflows/nf-core/bam_markduplicates_picard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info paramsSummaryLog(workflow)

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
// params.star_index       = getGenomeAttribute('star')
// params.rsem_index       = getGenomeAttribute('rsem')
params.hisat2_index     = getGenomeAttribute('hisat2')
// params.salmon_index     = getGenomeAttribute('salmon')
// params.kallisto_index   = getGenomeAttribute('kallisto')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description) : file("$projectDir/methods_description_template.yml", checkIfExists: true)

def multiqc_report     = []
def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

workflow {
ch_versions = Channel.empty()

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
ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)

// Prepare Genome
// TODO fasta empty
PREPARE_GENOME(
    params.fasta,
    params.gtf,
    params.splicesites,
    params.hisat2_index,
    params.gencode
)

// Aligning
ch_hisat2_multiqc = Channel.empty()
FASTQ_ALIGN_HISAT2(
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads, // reads       // channel: [ val(meta), [ reads ] ]
    PREPARE_GENOME.out.hisat2_index.map { [ [:], it ] }, // index       // channel: /path/to/hisat2/index
    PREPARE_GENOME.out.splicesites.map { [ [:], it ] }, // splicesites // channel: /path/to/genome.splicesites.txt
    PREPARE_GENOME.out.fasta.map { [ [:], it ] } // ch_fasta    // channel: [ fasta ]
)
ch_genome_bam        = FASTQ_ALIGN_HISAT2.out.bam
ch_genome_bam_index  = FASTQ_ALIGN_HISAT2.out.bai
ch_samtools_stats    = FASTQ_ALIGN_HISAT2.out.stats
ch_samtools_flagstat = FASTQ_ALIGN_HISAT2.out.flagstat
ch_samtools_idxstats = FASTQ_ALIGN_HISAT2.out.idxstats
ch_hisat2_multiqc    = FASTQ_ALIGN_HISAT2.out.summary
ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

// Markduplicates (optional)
ch_markduplicates_multiqc = Channel.empty()

BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
            PREPARE_GENOME.out.fai.map { [ [:], it ] }
        )
        ch_genome_bam             = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index       = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_samtools_stats         = BAM_MARKDUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = BAM_MARKDUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = BAM_MARKDUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = BAM_MARKDUPLICATES_PICARD.out.metrics
        //if (params.bam_csi_index) {
        //    ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.csi
        //}
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
//
// Get list of samples that failed trimming threshold for MultiQC report
//
ch_trim_read_count
    .map {
        meta, num_reads ->
            pass_trimmed_reads[meta.id] = true
            if (num_reads <= params.min_trimmed_reads.toFloat()) {
                pass_trimmed_reads[meta.id] = false
                return [ "$meta.id\t$num_reads" ]
            }
    }
    .collect()
    .map {
        tsv_data ->
            def header = ["Sample", "Reads after trimming"]
            multiqcTsvFromList(tsv_data, header)
    }
    .set { ch_fail_trimming_multiqc }

//
// MODULE: Pipeline reporting
//
CUSTOM_DUMPSOFTWAREVERSIONS (
    ch_versions.unique().collectFile(name: 'collated_versions.yml')
)

//
// MultiQC report (optional)
//
workflow_summary    = paramsSummaryMultiqc(summary_params)
ch_workflow_summary = Channel.value(workflow_summary)

methods_description    = methodsDescriptionText(ch_multiqc_custom_methods_description)
ch_methods_description = Channel.value(methods_description)

MULTIQC (
    ch_multiqc_config,
    ch_multiqc_custom_config.collect().ifEmpty([]),
    CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'),
    ch_multiqc_logo.collect().ifEmpty([]),
    ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]),
    //ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv').ifEmpty([]),
    //ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv').ifEmpty([]),
    ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
    ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
    ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_star_multiqc.collect{it[1]}.ifEmpty([]),
    ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_pseudo_multiqc.collect{it[1]}.ifEmpty([]),
    ch_samtools_stats.collect{it[1]}.ifEmpty([]),
    ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
    ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
    ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([])
    //ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_aligner_pca_multiqc.collect().ifEmpty([]),
    //ch_aligner_clustering_multiqc.collect().ifEmpty([]),
    //ch_pseudoaligner_pca_multiqc.collect().ifEmpty([]),
    //ch_pseudoaligner_clustering_multiqc.collect().ifEmpty([]),
    //ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_readduplication_multiqc.collect{it[1]}.ifEmpty([]),
    //ch_tin_multiqc.collect{it[1]}.ifEmpty([])
)
multiqc_report = MULTIQC.out.report.toList()
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