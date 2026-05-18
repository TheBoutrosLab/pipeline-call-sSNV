include { call_sSNV_Strelka2; call_sIndel_Manta } from './strelka2-processes'
include { filter_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common'

include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'

workflow strelka2 {
    take:
    META
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        call_sIndel_Manta(
            META,
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            params.reference_index,
            params.intersect_regions,
            params.intersect_regions_index
        )
        call_sSNV_Strelka2(
            META,
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            params.reference_index,
            call_sIndel_Manta.out[0],
            params.intersect_regions,
            params.intersect_regions_index
        )
        filter_VCF_BCFtools(META, call_sSNV_Strelka2.out.snvs_gzvcf
            .mix(call_sSNV_Strelka2.out.indels_gzvcf))
//  combine ids with each of the filtered strelka outputs (SNV and INDEL)
        Channel.from([['TUMOR', params.tumor_id], ['NORMAL', params.normal_id]])
            .map{ it -> ['orig_id': it[0], 'id': it[1]] }
            .collect()
            .combine(filter_VCF_BCFtools.out.gzvcf)
            .map { it.take(it.size() -2) }
            .set { rename_ids }
        rename_samples_BCFtools(
            META,
            rename_ids,
            filter_VCF_BCFtools.out.gzvcf
            )
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
        file_for_sha512 = indexed_vcfs
            .map{ it -> ["strelka2-${it[0]}-vcf", it[1]] }
            .mix( indexed_vcfs
            .map{ it -> ["strelka2-${it[0]}-index", it[2]] } )
        generate_sha512sum(META, file_for_sha512)
    emit:
        gzvcf = indexed_vcfs
            .filter { it[0] == 'SNV' }
            .map{ it -> ["${it[1]}"] }
        idx = indexed_vcfs
            .filter { it[0] == 'SNV' }
            .map{ it -> ["${it[2]}"] }
    }
