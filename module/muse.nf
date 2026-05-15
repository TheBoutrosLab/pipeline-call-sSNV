include { call_sSNV_MuSE; run_sump_MuSE } from './muse-processes'
include { filter_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'
workflow muse {
    take:
    META
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        call_sSNV_MuSE(
            META,
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            params.reference_index
        )
        run_sump_MuSE(
            META,
            call_sSNV_MuSE.out.txt,
            params.dbSNP,
            "${params.dbSNP}.tbi"
        )
        filter_VCF_BCFtools(META, run_sump_MuSE.out.vcf.map { it -> ['SNV', it] } )
        Channel.from([['TUMOR', params.tumor_id], ['NORMAL', params.normal_id]])
            .map{ it -> ['orig_id': it[0], 'id': it[1]] }
            .collect()
            .set { rename_ids }
        rename_samples_BCFtools(META, rename_ids, filter_VCF_BCFtools.out.gzvcf)
        compress_index_VCF(
            META.combine(rename_samples_BCFtools.out.gzvcf)
                .map{ it -> [
                    it[0] + [
                        "output_dir": it[0].workflow_output_dir,
                        "log_output_dir": "${it[0].log_output_dir}/process-log/${it[0].log_dir_prefix}",
                        "id": it[1],
                        "variant_type": it[1]
                    ],
                    it[2]
                ] }
            )
        indexed_vcfs = compress_index_VCF.out.index_out
            .map{ it -> [it[0].variant_type, it[1], it[2]] }
        file_for_sha512 = indexed_vcfs.map{ it -> ["muse-${it[0]}-vcf", it[1]] }
            .mix(indexed_vcfs.map{ it -> ["muse-${it[0]}-index", it[2]] })
        generate_sha512sum(META, file_for_sha512)
    emit:
        gzvcf = indexed_vcfs.map{ it -> ["${it[1]}"] }
        idx = indexed_vcfs.map{ it -> ["${it[2]}"] }
    }
