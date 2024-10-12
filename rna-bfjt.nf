// Import methods
include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet            } from 'plugin/nf-validation'
include { multiqcTsvFromList                                             } from './subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { paramsSummaryMultiqc                                           } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { getGenomeAttribute; multiqcTsvFromList; methodsDescriptionText } from './subworkflows/local/utils_nfcore_methods'

// Import modules
include { PICARD_MARKDUPLICATES                                          } from './modules/nf-core/picard/markduplicates/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                    } from './modules/nf-core/custom/dumpsoftwareversions/main'    
include { MULTIQC                                                        } from './modules/local/multiqc/main'                                 

// Import subworkflows
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE                               } from './subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { PREPARE_GENOME                                                 } from './subworkflows/local/prepare_genome/main'
include { FASTQ_ALIGN_HISAT2                                             } from './subworkflows/nf-core/fastq_align_hisat2/main'
include { BAM_MARKDUPLICATES_PICARD                                      } from './subworkflows/nf-core/bam_markduplicates_picard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def summary_params = paramsSummaryMap(workflow)
log.info paramsSummaryLog(workflow)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta            = getGenomeAttribute('fasta')
params.gff              = getGenomeAttribute('gff')
params.gtf              = getGenomeAttribute('gtf')
params.hisat2_index     = getGenomeAttribute('hisat2')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_multiqc_config                     = Channel.fromPath("$projectDir/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description) : file("$projectDir/methods_description_template.yml", checkIfExists: true)

def multiqc_report     = []
def pass_trimmed_reads = [:]

workflow {
ch_versions = Channel.empty()

//
// Reading the samplesheet file
//
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

//
// QC, UMI extraction, Adapter trimming
//
FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
    reads_ch,                                  // reads channel:                        [ val(meta), [ reads ] ]
    params.skip_fastqc,                        // skip fastqc:                          true/false
    params.with_umi,                           // enable UMI processing:                true/false
    params.skip_umi_extract,                   // skip UMI extraction:                  true/false
    params.skip_trimming,                      // skip adapter trimming:                true/false
    params.umi_discard_read,                   // read to discard post UMI extraction:  0, 1 or 2
    params.min_trimmed_reads                   // discard if trimmed reads fewer:       > 0
)
ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)

//
// Prepare Genome
//
PREPARE_GENOME(
    params.fasta,                              // fasta source file:           /path/to/genome.fasta
    params.gtf,                                // gtf source file:             /path/to/genome.gtf
    params.gff,                                // gff source file:             /path/to/genome.gff
    params.splicesites,                        // splicesites source file:     /path/to/splicesites.txt
    params.hisat2_index                        // hisat2 index directory:      /path/to/hisat2/index/
)
ch_hisat2_index = PREPARE_GENOME.out.hisat2_index
ch_splicesites = PREPARE_GENOME.out.splicesites
ch_fasta = PREPARE_GENOME.out.fasta
ch_fai = PREPARE_GENOME.out.fai
ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

//
// Alignment
//
ch_hisat2_multiqc = Channel.empty()

FASTQ_ALIGN_HISAT2(
    ch_filtered_reads,                         // reads channel:           [ val(meta), [ reads ] ]
    ch_hisat2_index.map { [ [:], it ] },       // index channel:           /path/to/hisat2/index
    ch_splicesites.map { [ [:], it ] },        // splicesites channel:     /path/to/genome.splicesites.txt
    ch_fasta.map { [ [:], it ] }               // fasta channel:           [ file(fasta) ]
)
ch_genome_bam        = FASTQ_ALIGN_HISAT2.out.bam
ch_genome_bam_index  = FASTQ_ALIGN_HISAT2.out.bai
ch_samtools_stats    = FASTQ_ALIGN_HISAT2.out.stats
ch_samtools_flagstat = FASTQ_ALIGN_HISAT2.out.flagstat
ch_samtools_idxstats = FASTQ_ALIGN_HISAT2.out.idxstats
ch_hisat2_multiqc    = FASTQ_ALIGN_HISAT2.out.summary
ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

//
// Markduplicates (if UMI extraction was skipped)
//
if (!params.skip_markduplicates && !params.with_umi) {
    ch_markduplicates_multiqc = Channel.empty()

    BAM_MARKDUPLICATES_PICARD (
                ch_genome_bam,                  // bam channel:     [ val(meta), path(reads) ]
                ch_fasta.map { [ [:], it ] },   // fasta channel:   [ path(fasta) ]
                ch_fai.map { [ [:], it ] }      // fai channel:     [ path(fai) ]
    )
    ch_genome_bam             = BAM_MARKDUPLICATES_PICARD.out.bam
    ch_genome_bam_index       = BAM_MARKDUPLICATES_PICARD.out.bai
    ch_samtools_stats         = BAM_MARKDUPLICATES_PICARD.out.stats
    ch_samtools_flagstat      = BAM_MARKDUPLICATES_PICARD.out.flagstat
    ch_samtools_idxstats      = BAM_MARKDUPLICATES_PICARD.out.idxstats
    ch_markduplicates_multiqc = BAM_MARKDUPLICATES_PICARD.out.metrics
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
}

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
// Get tool versions for MultiQC report
//
CUSTOM_DUMPSOFTWAREVERSIONS (
    ch_versions.unique().collectFile(name: 'collated_versions.yml')
)
ch_dumpversions = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml

//
// MultiQC report
//
workflow_summary    = paramsSummaryMultiqc(summary_params)
ch_workflow_summary = Channel.value(workflow_summary)

methods_description    = methodsDescriptionText(ch_multiqc_custom_methods_description)
ch_methods_description = Channel.value(methods_description)

MULTIQC (
    ch_multiqc_config,                                                                     
    ch_multiqc_custom_config.collect().ifEmpty([]),                                        
    ch_dumpversions.collect(),                                                             
    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),                    
    ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'),              
    ch_multiqc_logo.collect().ifEmpty([]),                                                 
    ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]),
    ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),                                      
    ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),                                     
    ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),                                        
    ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),                                          
    ch_samtools_stats.collect{it[1]}.ifEmpty([]),                                          
    ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),                                       
    ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),                                       
    ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([])                                   
)
multiqc_report = MULTIQC.out.report.toList()
}
