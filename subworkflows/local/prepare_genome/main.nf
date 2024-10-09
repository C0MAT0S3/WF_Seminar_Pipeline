//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SORTMERNA_INDEX    } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_KALLISTO_INDEX     } from '../../../modules/nf-core/untar'

include { HISAT2_EXTRACTSPLICESITES         } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                      } from '../../../modules/nf-core/hisat2/build'

include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { GTF_FILTER                           } from '../../../modules/local/gtf_filter'

workflow PREPARE_GENOME {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    //gff                      //      file: /path/to/genome.gff
    //additional_fasta         //      file: /path/to/additional.fasta
    //transcript_fasta         //      file: /path/to/transcript.fasta
    //gene_bed                 //      file: /path/to/gene.bed
    splicesites              //      file: /path/to/splicesites.txt
    //bbsplit_fasta_list       //      file: /path/to/bbsplit_fasta_list.txt
    //sortmerna_fasta_list     //      file: /path/to/sortmerna_fasta_list.txt
    //star_index               // directory: /path/to/star/index/
    //rsem_index               // directory: /path/to/rsem/index/
    //salmon_index             // directory: /path/to/salmon/index/
    //kallisto_index           // directory: /path/to/kallisto/index/
    hisat2_index             // directory: /path/to/hisat2/index/
    //bbsplit_index            // directory: /path/to/bbsplit/index/
    //sortmerna_index          // directory: /path/to/sortmerna/index/
    gencode                  //   boolean: whether the genome is from GENCODE
    //featurecounts_group_type //    string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    //aligner                  //    string: Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem' and 'hisat2'
    //pseudo_aligner           //    string: Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'
    //skip_gtf_filter          //   boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    //skip_bbsplit             //   boolean: Skip BBSplit for removal of non-reference genome reads
    //skip_sortmerna           //   boolean: Skip sortmerna for removal of reads mapping to sequences in sortmerna_fasta_list
    //skip_alignment           //   boolean: Skip all of the alignment-based processes within the pipeline
    //skip_pseudo_alignment    //   boolean: Skip all of the pseudoalignment-based processes within the pipeline

    main:
    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], file(fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
    }

     //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], file(gtf, checkIfExists: true) ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value(file(gtf, checkIfExists: true))
            }
        } else if (gff) {
            if (gff.endsWith('.gz')) {
                ch_gff      = GUNZIP_GFF ( [ [:], file(gff, checkIfExists: true) ] ).gunzip
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = Channel.value(file(gff, checkIfExists: true)).map { [ [:], it ] }
            }
            ch_gtf      = GFFREAD ( ch_gff, [] ).gtf.map { it[1] }
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if (!splicesites) {
        ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
        ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
    } else {
        ch_splicesites = Channel.value(file(splicesites))
    }
    if (hisat2_index) {
        if (hisat2_index.endsWith('.tar.gz')) {
            ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], hisat2_index ] ).untar.map { it[1] }
            ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
        } else {
            ch_hisat2_index = Channel.value(file(hisat2_index))
        }
    } else {
        ch_hisat2_index = HISAT2_BUILD ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
        ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
    }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    //fai              = ch_fai                    // channel: path(genome.fai)
    //gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    //transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    //chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    //bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    //rrna_fastas      = ch_rrna_fastas            // channel: path(sortmerna_fasta_list)
    //sortmerna_index  = ch_sortmerna_index        // channel: path(sortmerna/index/)
    //star_index       = ch_star_index             // channel: path(star/index/)
    //rsem_index       = ch_rsem_index             // channel: path(rsem/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    //salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    //kallisto_index   = ch_kallisto_index         // channel: [ meta, path(kallisto/index/) ]
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}