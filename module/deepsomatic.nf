include { run_SplitIntervals_GATK; run_MergeVcfs_GATK } from './mutect2-processes'
include { convert_IntervalListToBed_GATK; call_sSNV_DeepSomatic } from './deepsomatic-processes'
include { filter_VCF_BCFtools; split_VCF_BCFtools; generate_sha512sum } from './common' addParams(
    log_dir_prefix: "DeepSomatic-${params.deepsomatic_version}"
    )
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: "${params.log_output_dir}/process-log/DeepSomatic-${params.deepsomatic_version}",
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
        ])

workflow deepsomatic {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        run_SplitIntervals_GATK(
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
                    file(interval_path).getName().replace('-contig.interval_list', ''),
                    interval_path
                ]
            } | convert_IntervalListToBed_GATK


        call_sSNV_DeepSomatic(
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
            .collect() | run_MergeVcfs_GATK

        run_MergeVcfs_GATK.out.unfiltered
            .map{ unfiltered_vcf -> ['all', unfiltered_vcf] } | filter_VCF_BCFtools

        split_VCF_BCFtools(
            filter_VCF_BCFtools.out.gzvcf.map{ filtered_vcf -> filtered_vcf[1] },
            ['snps', 'mnps', 'indels']
        )

        compress_index_VCF(
            split_VCF_BCFtools.out.gzvcf
        )

        compress_index_VCF.out.index_out
            .map{ indexed_out ->
                ["deepsomatic-${indexed_out[0]}-vcf", indexed_out[1]]
            }
            .mix(
                compress_index_VCF.out.index_out
                    .map{ index_file ->
                        ["deepsomatic-${index_file[0]}-index", index_file[2]]
                    }
            )
            .set{ files_for_checksum }

        generate_sha512sum(files_for_checksum)
    }
