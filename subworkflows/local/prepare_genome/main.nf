//
// Uncompress and prepare reference genome files
//
include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../../modules/nf-core/untar'

include { HISAT2_EXTRACTSPLICESITES         } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                      } from '../../../modules/nf-core/hisat2/build'

include { CUSTOM_GETCHROMSIZES              } from '../../../modules/nf-core/custom/getchromsizes'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'

workflow PREPARE_GENOME {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    gff                      //      file: /path/to/genome.gff
    splicesites              //      file: /path/to/splicesites.txt
    hisat2_index             // directory: /path/to/hisat2/index/

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
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

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
    fai              = ch_fai                    // channel: path(genome.fai)
    splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}