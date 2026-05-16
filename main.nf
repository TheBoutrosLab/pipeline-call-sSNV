#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference_index = "${params.reference}.fai"
params.reference_dict = "${file(params.reference).parent / file(params.reference).baseName}.dict"

include { generate_standard_filename } from './external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include { run_validate_PipeVal } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf'
include { indexFile } from './external/pipeline-Nextflow-module/modules/common/indexFile/main.nf'
include { somaticsniper } from './module/somaticsniper'
include { strelka2 } from './module/strelka2'
include { mutect2 } from './module/mutect2'
include { muse } from './module/muse'
include { deepsomatic } from './module/deepsomatic'
include { process_vcfs } from './module/process-vcfs'
include { intersect; getToolName } from './module/intersect'
include { plot_vaf } from './module/plot-vaf'

log.info """\
    ------------------------------------
    C A L L - S S N V    P I P E L I N E
    ------------------------------------
    Boutros Lab

    Current Configuration:
    - pipeline:
        name: ${workflow.manifest.name}
        version: ${workflow.manifest.version}

    - input:
        dataset_id: ${params.dataset_id}
        patient_id: ${params.patient_id}
        sample_id: ${params.sample_id}
        algorithm: ${params.algorithm}
        tumor: ${params.samples_to_process.findAll{ it.sample_type == 'tumor' }['path']}
        normal: ${params.samples_to_process.findAll{ it.sample_type == 'normal' }['path']}
        reference: ${params.reference}
        intersect_regions: ${params.intersect_regions}

    - output:
        output_dir: ${params.output_dir_base}
        log_output_dir: ${params.log_output_dir}

    - option:
        ucla_cds: ${params.ucla_cds}
        save_intermediate_files: ${params.save_intermediate_files}
        docker_container_registry: ${params.docker_container_registry}
        bgzip_extra_args = ${params.bgzip_extra_args}
        tabix_extra_args = ${params.tabix_extra_args}

    - sample names extracted from input BAM files and sanitized:
        tumor_in: ${params.samples_to_process.findAll{ it.sample_type == 'tumor' }['orig_id']}
        tumor_out: ${params.samples_to_process.findAll{ it.sample_type == 'tumor' }['id']}
        normal_in: ${params.samples_to_process.findAll{ it.sample_type == 'normal' }['orig_id']}
        normal_out: ${params.samples_to_process.findAll{ it.sample_type == 'normal' }['id']}
"""

if (params.input_type == 'bam') {
    if (params.max_cpus < 8 || params.max_memory < 16) {
        if (params.algorithm.contains('muse') || params.algorithm.contains('mutect2')) {
            throw new Exception(
                "Insufficient resources: ${params.max_cpus} CPUs and ${params.max_memory} of memory." +
                " To run Mutect2 this pipeline requires at least 8 CPUs and 16 GB of memory." +
                " To run MuSE this pipeline requires at least 16 CPUs and 32 GB of memory."
                )
            }
        }
    else if (params.max_cpus < 16 || params.max_memory < 32) {
        if (params.algorithm.contains('muse')) {
            throw new Exception(
                "Insufficient resources: ${params.max_cpus} CPUs and ${params.max_memory} of memory." +
                " To run MuSE this pipeline requires at least 16 CPUs and 32 GB of memory."
                )
            }
        }

    Channel
        .from( params.samples_to_process )
            .filter{ it.sample_type == 'tumor' }
            .multiMap{ it ->
                tumor_bam: it['path']
                tumor_index: indexFile(it['path'])
                contamination_est: it['contamination_table']
                }
            .set { tumor_input_chs }

    Channel
        .from( params.samples_to_process )
            .filter{ it.sample_type == 'normal' }
            .ifEmpty(['path': "${params.work_dir}/NO_FILE.bam"])
            .multiMap{ it ->
                normal_bam: it['path']
                normal_index: indexFile(it['path'])
                }
            .set { normal_input_chs }

    // Set empty channels so any unused tools don't cause failure at intersect step
    Channel.empty().set { somaticsniper_gzvcf_ch }
    Channel.empty().set { strelka2_gzvcf_ch }
    Channel.empty().set { muse_gzvcf_ch }
    Channel.empty().set { mutect2_gzvcf_ch }
    Channel.empty().set { deepsomatic_gzvcf_ch }

    Channel.empty().set { somaticsniper_idx_ch }
    Channel.empty().set { strelka2_idx_ch }
    Channel.empty().set { muse_idx_ch }
    Channel.empty().set { mutect2_idx_ch }
    Channel.empty().set { deepsomatic_idx_ch }

} else if (params.input_type == 'vcf') {
    Channel
        .fromList(params.samples_to_process)
        .map { vcf ->
            return tuple(vcf.path, indexFile(vcf.path), vcf.algorithm)
        }
        .set { samplesToProcess_ch }
    }

script_dir_ch = Channel.fromPath(
    "$projectDir/script",
    checkIfExists: true
    )

workflow {
    meta_base = Channel.value([
        "output_dir_base": params.output_dir_base,
        "log_output_dir": params.log_output_dir,
        "save_intermediate_files": params.save_intermediate_files,
        "bgzip_extra_args": params.bgzip_extra_args,
        "tabix_extra_args": params.tabix_extra_args
    ])

    pipeval_meta = meta_base.map{ base_m ->
        base_m + [
            "docker_image": params.docker_image_validate_params
        ]
    }

    somaticsniper_meta = meta_base.map{ base_m ->
        base_m + [
            "workflow_output_dir": "${params.output_dir_base}/SomaticSniper-${params.somaticsniper_version}",
            "output_filename": generate_standard_filename(
                "SomaticSniper-${params.somaticsniper_version}",
                params.dataset_id,
                params.sample_id,
                [:]
                ),
            "log_dir_prefix": "SomaticSniper-${params.somaticsniper_version}"
        ]
    }

    strelka2_meta = meta_base.map{ base_m ->
        base_m + [
            "workflow_output_dir": "${params.output_dir_base}/Strelka2-${params.strelka2_version}",
            "output_filename": generate_standard_filename(
                "Strelka2-${params.strelka2_version}",
                params.dataset_id,
                params.sample_id,
                [:]
                ),
            "log_dir_prefix": "Strelka2-${params.strelka2_version}"
        ]
    }

    mutect2_meta = meta_base.map{ base_m ->
        base_m + [
            "workflow_output_dir": "${params.output_dir_base}/Mutect2-${params.GATK_version}",
            "output_filename": generate_standard_filename(
                "Mutect2-${params.GATK_version}",
                params.dataset_id,
                params.sample_id,
                [:]
                ),
            "log_dir_prefix": "Mutect2-${params.GATK_version}"
        ]
    }

    muse_meta = meta_base.map{ base_m ->
        base_m + [
            "workflow_output_dir": "${params.output_dir_base}/MuSE-${params.MuSE_version}",
            "output_filename": generate_standard_filename(
                "MuSE-${params.MuSE_version}",
                params.dataset_id,
                params.sample_id,
                [:]
                ),
            "log_dir_prefix": "MuSE-${params.MuSE_version}"
        ]
    }

    deepsomatic_meta = meta_base.map{ base_m ->
        base_m + [
            "workflow_output_dir": "${params.output_dir_base}/DeepSomatic-${params.deepsomatic_version}",
            "output_filename": generate_standard_filename(
                "DeepSomatic-${params.deepsomatic_version}",
                params.dataset_id,
                params.sample_id,
                [:]
                ),
            "log_dir_prefix": "DeepSomatic-${params.deepsomatic_version}",
            "unzip_and_rezip": true
        ]
    }

    intersect_meta = meta_base.map{ base_m ->
        base_m + [
            "workflow_output_dir": "${params.output_dir_base}/Intersect-BCFtools-${params.BCFtools_version}",
            "output_filename": generate_standard_filename(
                "BCFtools-${params.BCFtools_version}",
                params.dataset_id,
                params.sample_id,
                [:]
                ),
            "log_dir_prefix": "Intersect-BCFtools-${params.BCFtools_version}",
            "compress_publishdir": "${params.output_dir_base}/Intersect-BCFtools-${params.BCFtools_version}/output",
            "compress_enabled": true
        ]
    }

    process_vcfs_meta = intersect_meta.map{ base_m ->
        base_m + [
            "keep_input_prefix": true
        ]
    }

    plot_vaf_meta = intersect_meta.map{ base_m ->
        base_m + [
            "output_filename": generate_standard_filename(
                "BPG-${params.bpg_version}",
                params.dataset_id,
                params.sample_id,
                [:]
                )
        ]
    }

    reference_ch = Channel.from(
        params.reference,
        params.reference_index,
        params.reference_dict
        )

    intersect_regions_ch = Channel.from(
        params.intersect_regions,
        params.intersect_regions_index
        )

    if (params.input_type == 'bam') {
        files_to_validate_ch = reference_ch
            .mix(intersect_regions_ch)
            .mix(tumor_input_chs.tumor_bam, tumor_input_chs.tumor_index)

        if (params.samples_to_process.findAll{ it.sample_type == 'normal' }.size() > 0) {
            files_to_validate_ch = files_to_validate_ch
                .mix(normal_input_chs.normal_bam, normal_input_chs.normal_index)
            }

        run_validate_PipeVal(pipeval_meta.combine(files_to_validate_ch))
        run_validate_PipeVal.out.validation_result.collectFile(
            name: 'input_validation.txt', newLine: true,
            storeDir: "${params.output_dir_base}/validation"
            )

        if ('somaticsniper' in params.algorithm) {
            somaticsniper(
                somaticsniper_meta,
                tumor_input_chs.tumor_bam,
                tumor_input_chs.tumor_index,
                normal_input_chs.normal_bam,
                normal_input_chs.normal_index
            )
            somaticsniper.out.gzvcf.set { somaticsniper_gzvcf_ch }
            somaticsniper.out.idx.set { somaticsniper_idx_ch }
            }
        if ('strelka2' in params.algorithm) {
            strelka2(
                strelka2_meta,
                tumor_input_chs.tumor_bam,
                tumor_input_chs.tumor_index,
                normal_input_chs.normal_bam,
                normal_input_chs.normal_index
            )
            strelka2.out.gzvcf.set { strelka2_gzvcf_ch }
            strelka2.out.idx.set { strelka2_idx_ch }
            }
        if ('muse' in params.algorithm) {
            muse(
                muse_meta,
                tumor_input_chs.tumor_bam,
                tumor_input_chs.tumor_index,
                normal_input_chs.normal_bam,
                normal_input_chs.normal_index
            )
            muse.out.gzvcf.set { muse_gzvcf_ch }
            muse.out.idx.set { muse_idx_ch }
            }
        if ('mutect2' in params.algorithm) {
            mutect2(
                mutect2_meta,
                tumor_input_chs.tumor_bam.collect(),
                tumor_input_chs.tumor_index.collect(),
                normal_input_chs.normal_bam.collect(),
                normal_input_chs.normal_index.collect(),
                tumor_input_chs.contamination_est.collect()
            )
            mutect2.out.gzvcf.set { mutect2_gzvcf_ch }
            mutect2.out.idx.set { mutect2_idx_ch }
            }
        if ('deepsomatic' in params.algorithm) {
            deepsomatic(
                deepsomatic_meta,
                tumor_input_chs.tumor_bam,
                tumor_input_chs.tumor_index,
                normal_input_chs.normal_bam,
                normal_input_chs.normal_index
            )
            deepsomatic.out.gzvcf.set{ deepsomatic_gzvcf_ch }
            deepsomatic.out.idx.set{ deepsomatic_idx_ch }
        }

        tool_gzvcfs = somaticsniper_gzvcf_ch
            .mix(strelka2_gzvcf_ch)
            .mix(mutect2_gzvcf_ch)
            .mix(muse_gzvcf_ch)
            .mix(deepsomatic_gzvcf_ch)
            .collect()
        tool_indices = somaticsniper_idx_ch
            .mix(strelka2_idx_ch)
            .mix(mutect2_idx_ch)
            .mix(muse_idx_ch)
            .mix(deepsomatic_idx_ch)
            .collect()

    } else if (params.input_type == 'vcf') {
        process_vcfs(process_vcfs_meta, samplesToProcess_ch)
        process_vcfs.out.gzvcf.set { tool_gzvcfs }
        process_vcfs.out.idx.set { tool_indices }
        }

    if (params.algorithm.size() > 1 || params.input_type == 'vcf') {
        intersect(
            intersect_meta,
            tool_gzvcfs,
            tool_indices,
            script_dir_ch,
            )
        }

    all_files = tool_gzvcfs.mix(tool_indices)
        .flatten()
        .collect()

    identified_gzvcfs = tool_gzvcfs.flatten()
        .map{ [algorithm: getToolName(it), path: it] }
        .collect()

    plot_vaf(
        plot_vaf_meta,
        identified_gzvcfs,
        all_files
        )

    }
