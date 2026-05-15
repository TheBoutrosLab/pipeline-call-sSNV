include { run_SplitIntervals_GATK; run_MergeVcfs_GATK } from './mutect2-processes'
include { convert_IntervalListToBed_GATK; call_sSNV_DeepSomatic } from './deepsomatic-processes'
include { filter_VCF_BCFtools; split_VCF_BCFtools; generate_sha512sum } from './common'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'

workflow deepsomatic {
    take:
    META
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        run_SplitIntervals_GATK(
            META,
            params.intersect_regions,
            params.intersect_regions_index,
            params.reference,
            params.reference_index,
            params.reference_dict
        )

        run_SplitIntervals_GATK.out.interval_list
            .flatten()
            .map{ interval_path ->
                [
                    file(interval_path).getName().replace('-scattered.interval_list', ''),
                    interval_path
                ]
            }
            .set{ interval_list_bed_input }

        convert_IntervalListToBed_GATK(
            META,
            interval_list_bed_input
        )


        call_sSNV_DeepSomatic(
            META,
            convert_IntervalListToBed_GATK.out.interval_bed,
            tumor_bam
                .combine(convert_IntervalListToBed_GATK.out.interval_bed)
                .map{ tumor_bam_combined -> tumor_bam_combined[0] },
            tumor_index
                .combine(convert_IntervalListToBed_GATK.out.interval_bed)
                .map{ tumor_index_combined -> tumor_index_combined[0] },
            normal_bam
                .combine(convert_IntervalListToBed_GATK.out.interval_bed)
                .map{ normal_bam_combined -> normal_bam_combined[0] },
            normal_index
                .combine(convert_IntervalListToBed_GATK.out.interval_bed)
                .map{ normal_index_combined -> normal_index_combined[0] },
            params.reference,
            params.reference_index,
            params.reference_dict,
        )

        call_sSNV_DeepSomatic.out.vcf
            .filter{ raw_vcf_out -> (raw_vcf_out[2] != "0") } // Filter out empty VCFs
            .map{ filtered_vcf -> filtered_vcf[0] }
            .collect()
            .set{ deepsomatic_vcfs_for_merge }

        run_MergeVcfs_GATK(
            META,
            deepsomatic_vcfs_for_merge
        )

        run_MergeVcfs_GATK.out.unfiltered
            .map{ unfiltered_vcf -> ['all', unfiltered_vcf] }
            .set{ filtered_vcf_input }

        filter_VCF_BCFtools(
            META,
            filtered_vcf_input
        )

        split_VCF_BCFtools(
            META,
            filter_VCF_BCFtools.out.gzvcf.map{ filtered_vcf -> filtered_vcf[1] },
            ['snps', 'mnps', 'indels']
        )

        compress_index_VCF(
            META.combine(split_VCF_BCFtools.out.gzvcf)
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

        indexed_vcfs
            .map{ indexed_out ->
                ["deepsomatic-${indexed_out[0]}-vcf", indexed_out[1]]
            }
            .mix(
                indexed_vcfs
                    .map{ index_file ->
                        ["deepsomatic-${index_file[0]}-index", index_file[2]]
                    }
            )
            .set{ files_for_checksum }

        generate_sha512sum(META, files_for_checksum)

    emit:
        gzvcf = indexed_vcfs
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[1]}"] }
        idx = indexed_vcfs
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[2]}"] }
}
