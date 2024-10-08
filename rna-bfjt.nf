// Local modules
include { NCBIGENOMEDOWNLOAD } '../modules/nf-core/ncbigenomedownload/main'
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main' 
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { MARK_DUPLICATES } from '../modules/nf-core/picard/markduplicates/main'

workflow {
// Reading
TODO
// FastQC
TODO
// Trimming
TODO
// Aligning
TODO
// MultiQC report (optional)
TODO
// Markduplicates (optional)
TODO
}