include { run_SplitIntervals_GATK; call_sSNV_Mutect2; run_MergeVcfs_GATK; run_MergeMutectStats_GATK; run_LearnReadOrientationModel_GATK; run_FilterMutectCalls_GATK } from './mutect2-processes'
include { filter_VCF_BCFtools; split_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'

workflow mutect2 {
    take:
    META
    tumor_bam
    tumor_index
    normal_bam
    normal_index
    contamination_table

    main:
        run_SplitIntervals_GATK(
            META,
            params.intersect_regions,
            params.intersect_regions_index,
            params.reference,
            params.reference_index,
            params.reference_dict
            )

        Channel
            .from( params.samples_to_process )
            .map{ it -> ['orig_id': it['orig_id'], 'id': it['id'], 'sample_type': it['sample_type']] }
            .set { id_ch }

        normal_orig_ids = id_ch
            .filter{ it['sample_type'] == 'normal' }
            .ifEmpty(['orig_id': 'NO_ID'])
            .map{ it['orig_id'] }
            .collect()

        // to avoid input file name collision or null input error in Mutect2
        contamination_table
            .flatten()
            .unique()
            .filter{ it !== null }
            .ifEmpty("${params.work_dir}/NO_PATH")
            .set { contamination_table }

        call_sSNV_Mutect2(
            META,
            run_SplitIntervals_GATK.out.interval_list.flatten(),
            tumor_bam
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            tumor_index
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            normal_bam
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            normal_index
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            params.reference,
            params.reference_index,
            params.reference_dict,
            normal_orig_ids
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            params.germline_resource_gnomad_vcf,
            params.germline_resource_gnomad_vcf_index
            )

        ich_MergeVcfs = call_sSNV_Mutect2.out.unfiltered.collect()
        ich_MergeMutectStats = call_sSNV_Mutect2.out.unfiltered_stats.collect()
        ich_LearnReadOrientationModel = call_sSNV_Mutect2.out.f1r2.collect()

        run_MergeVcfs_GATK(META, ich_MergeVcfs)
        run_MergeMutectStats_GATK(META, ich_MergeMutectStats)
        run_LearnReadOrientationModel_GATK(META, ich_LearnReadOrientationModel)
        run_FilterMutectCalls_GATK(
            META,
            params.reference,
            params.reference_index,
            params.reference_dict,
            run_MergeVcfs_GATK.out.unfiltered,
            run_MergeVcfs_GATK.out.unfiltered_index,
            run_MergeMutectStats_GATK.out.merged_stats,
            run_LearnReadOrientationModel_GATK.out.read_orientation_model,
            contamination_table.collect()
            )
        filter_VCF_BCFtools(META, run_FilterMutectCalls_GATK.out.filtered
            .map{ it -> ['all', it] }
            )
        split_VCF_BCFtools(META, filter_VCF_BCFtools.out.gzvcf
            .map{ it -> it[1] },
            ['snps', 'mnps', 'indels']
        )
        rename_samples_BCFtools(
            META,
            // combine with split_VCF_BCFtools output to duplicate the id input for each file.
            id_ch
                .collect()
                .combine(split_VCF_BCFtools.out.gzvcf)
                .map { it.take(it.size() -2) } //remove the split_VCF_BCFtools files
            ,
            split_VCF_BCFtools.out.gzvcf
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
            .map{ it -> ["mutect2-${it[0]}-vcf", it[1]] }
            .mix( indexed_vcfs
            .map{ it -> ["mutect2-${it[0]}-index", it[2]] }
            )
        generate_sha512sum(META, file_for_sha512)
    emit:
        gzvcf = indexed_vcfs
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[1]}"] }
        idx = indexed_vcfs
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[2]}"] }
    }
