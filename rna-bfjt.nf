// Import modules
include { NCBIGENOMEDOWNLOAD } from './modules/nf-core/ncbigenomedownload/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'
include { STAR_ALIGN } from './modules/nf-core/star/align/main' 
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main'

workflow {
// Reading the samplesheet file
samplesheet_file = file(params.input)
samplesheet_dir = samplesheet_file.parent

// Read the samplesheet into a channel
input_ch = channel.fromPath(params.input).splitCsv(header: true)

// Create meta map for analyses
    map_ch = input_ch.map { row ->
        // Resolve relative fastq files in the samplesheet
        def fastq1 = file("${samplesheet_dir}/${row.fastq_1}")
        def fastq2 = file("${samplesheet_dir}/${row.fastq_2}")
        [row.subMap("sample", "strandedness"), [fastq1, fastq2]]
    }
map_ch.groupTuple()

// FastQC
FASTQC(map_ch)
// Trimming
//TODO
// Aligning
//TODO
// MultiQC report (optional)
//TODO
// Markduplicates (optional)
//TODO
}