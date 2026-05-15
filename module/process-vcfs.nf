include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include { split_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'

workflow process_vcfs {
    take:
        META
        samplesToProcess_ch
    main:
        split_VCF_BCFtools(META, samplesToProcess_ch
            .filter { it[2] == 'mutect2' }
            .map{ it -> it[0] },
            ['snps', 'mnps', 'indels']
            )
        rename_files_ch = samplesToProcess_ch
            .filter { it[2] != 'mutect2' }
            .map{ it -> ['SNV', it[0]] }
            .mix(split_VCF_BCFtools.out.gzvcf)
        rename_id_ch = Channel.value(['orig_id': params.input_tumor_id,'id': params.tumor_id, 'sample_type': 'tumor' ])
            .mix(Channel.value(['orig_id': params.input_normal_id, 'id': params.normal_id, 'sample_type': 'normal' ]))
            .mix(Channel.value(['orig_id': 'TUMOR', 'id': params.tumor_id, 'sample_type': 'tumor' ]))
            .mix(Channel.value(['orig_id': 'NORMAL', 'id': params.normal_id, 'sample_type': 'normal' ]))
            .collect()
        rename_samples_BCFtools(
            META,
            rename_id_ch
                .combine(rename_files_ch) // combine with files channel to get the count right
                .map { it.take(it.size() -2) } //remove the files
            ,
            rename_files_ch
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
            .map{ it -> ["${file(it[0]).getName().split('-')[0]}-vcf", it[1]] }
            .mix( indexed_vcfs
            .map{ it -> ["${file(it[0]).getName().split('-')[0]}-index", it[2]] }
            )
        generate_sha512sum(META, file_for_sha512)

    emit:
        gzvcf = indexed_vcfs
            .filter { it[0].replace('snps', 'SNV') == 'SNV' }
            .map{ it -> ["${it[1]}"] }
        idx = indexed_vcfs
            .filter { it[0].replace('snps', 'SNV') == 'SNV' }
            .map{ it -> ["${it[2]}"] }
    }
