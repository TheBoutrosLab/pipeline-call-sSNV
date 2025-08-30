include { run_SplitIntervals_GATK; run_MergeVcfs_GATK } from './mutect2-processes'
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
                'interval_id': file(interval_path).getName().replace('-contig.interval_list', ''),
                'interval_path': interval_path
            ]
        }
        .set{ input_ch_intervals }

        filter_VCF_BCFtools(run_FilterMutectCalls_GATK.out.filtered
            .map{ it -> ['all', it] }
            )

        split_VCF_BCFtools(filter_VCF_BCFtools.out.gzvcf
            .map{ it -> it[1] },
            ['snps', 'mnps', 'indels']
        )

        compress_index_VCF(rename_samples_BCFtools.out.gzvcf)

        file_for_sha512 = compress_index_VCF.out.index_out
            .map{ it -> ["deepsomatic-${it[0]}-vcf", it[1]] }
            .mix( compress_index_VCF.out.index_out
            .map{ it -> ["deepsomatic-${it[0]}-index", it[2]] }
            )
        generate_sha512sum(file_for_sha512)

    emit:
        gzvcf = compress_index_VCF.out.index_out
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[1]}"] }
        idx = compress_index_VCF.out.index_out
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[2]}"] }
    }
