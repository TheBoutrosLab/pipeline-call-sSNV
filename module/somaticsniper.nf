include { call_sSNV_SomaticSniper; convert_BAM2Pileup_SAMtools; create_IndelCandidate_SAMtools; apply_NormalIndelFilter_SomaticSniper; apply_TumorIndelFilter_SomaticSniper; create_ReadCountPosition_SomaticSniper; generate_ReadCount_bam_readcount; filter_FalsePositive_SomaticSniper; call_HighConfidenceSNV_SomaticSniper } from './somaticsniper-processes'
include { rename_samples_BCFtools; generate_sha512sum } from './common'
include { compress_index_VCF as compress_index_VCF_hc } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'
include { compress_index_VCF as compress_index_VCF_fix } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'
include { compress_file_bzip2 } from './common'

workflow somaticsniper {
    take:
    META
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        call_sSNV_SomaticSniper(
            META,
            tumor_bam,
            normal_bam,
            params.reference,
            params.reference_index
            )
        tumor_bam_path = tumor_bam
            .map{it -> ['tumor', it]}
        normal_bam_path = normal_bam
            .map{it -> ['normal', it]}
        ch_convert_BAM2Pileup_SAMtools_bams = tumor_bam_path
            .mix(normal_bam_path)
        convert_BAM2Pileup_SAMtools(
            META,
            ch_convert_BAM2Pileup_SAMtools_bams,
            params.reference,
            params.reference_index
            )
        create_IndelCandidate_SAMtools(META, convert_BAM2Pileup_SAMtools.out.raw_pileup)

        // tumor and normal need to be processed seperately.
        create_IndelCandidate_SAMtools.out.filtered_pileup
            .branch {
                normal: it[0] == "normal"
                        return it[1]
                tumor: it[0] == "tumor"
                        return it[1]
            }
            .set { ch_snpfilter }

        apply_NormalIndelFilter_SomaticSniper(
            META,
            call_sSNV_SomaticSniper.out.bam_somaticsniper,
            ch_snpfilter.normal
            )
        apply_TumorIndelFilter_SomaticSniper(
            META,
            apply_NormalIndelFilter_SomaticSniper.out.vcf_normal,
            ch_snpfilter.tumor
            )
        create_ReadCountPosition_SomaticSniper(
            META,
            apply_TumorIndelFilter_SomaticSniper.out.vcf_tumor
            )
        generate_ReadCount_bam_readcount(
            META,
            params.reference,
            params.reference_index,
            create_ReadCountPosition_SomaticSniper.out.snp_positions,
            tumor_bam,tumor_index
            )
        filter_FalsePositive_SomaticSniper(
            META,
            apply_TumorIndelFilter_SomaticSniper.out.vcf_tumor,
            generate_ReadCount_bam_readcount.out.readcount
            )
        call_HighConfidenceSNV_SomaticSniper(
            META,
            filter_FalsePositive_SomaticSniper.out.fp_pass
            )
        // combining to delay compression until after filtering step
        compress_file_bzip2(
            META.map{ base_m ->
                base_m + [
                    "compress_publishdir": "${base_m.workflow_output_dir}/intermediate/generate_ReadCount_bam_readcount",
                    "compress_enabled": base_m.save_intermediate_files
                ]
            },
            generate_ReadCount_bam_readcount.out.readcount
                .combine(filter_FalsePositive_SomaticSniper.out.fp_pass.collect())
                .map{ it -> ['readcount', it[0]] }
            )
        // rename_samples_BCFtools needs bgzipped input
        compress_index_VCF_hc(
            META.combine(call_HighConfidenceSNV_SomaticSniper.out.hc_vcf.map{ it -> ['SNV', it] })
                .map{ it -> [
                    it[0] + [
                        "output_dir": it[0].workflow_output_dir,
                        "log_output_dir": "${it[0].log_output_dir}/process-log/${it[0].log_dir_prefix}",
                        "id": it[1],
                        "variant_type": it[1],
                        "is_output_file": false
                    ],
                    it[2]
                ] }
            )
        hc_indexed_vcfs = compress_index_VCF_hc.out.index_out
            .map{ it -> [it[0].variant_type, it[1], it[2]] }
        Channel.from([['TUMOR', params.tumor_id], ['NORMAL', params.normal_id]])
            .map{ it -> ['orig_id': it[0], 'id': it[1]] }
            .collect()
            .set { rename_ids }
        rename_samples_BCFtools(META, rename_ids, hc_indexed_vcfs
            .map{ it -> [it[0], it[1]] })
        compress_index_VCF_fix(
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
        fixed_indexed_vcfs = compress_index_VCF_fix.out.index_out
            .map{ it -> [it[0].variant_type, it[1], it[2]] }
        file_for_sha512 = fixed_indexed_vcfs
            .map{ it -> ["${it[0]}-vcf", it[1]] }
            .mix(fixed_indexed_vcfs
                .map{ it -> ["${it[0]}-index", it[2]] }
                )
        generate_sha512sum(META, file_for_sha512)
    emit:
        gzvcf = fixed_indexed_vcfs.map{ it -> ["${it[1]}"] }
        idx = fixed_indexed_vcfs.map{ it -> ["${it[2]}"] }
    }
