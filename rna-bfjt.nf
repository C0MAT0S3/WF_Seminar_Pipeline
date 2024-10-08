// Local modules
include { NCBIGENOMEDOWNLOAD } from './modules/nf-core/ncbigenomedownload/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { STAR_ALIGN } from './modules/nf-core/star/align/main' 
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main'

workflow {
// Reading
input_ch = channel.fromPath(params.input).splitCsv(header: true)

// Create meta map for analyses
map_ch = input_ch.map { row -> [row.subMap("sample", "strandedness"), [file(row.fastq_1), file(row.fastq_2)]] }
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