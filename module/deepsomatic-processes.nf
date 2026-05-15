log.info """\
====================================
        D E E P S O M A T I C
====================================
Docker Images:
- docker_image_GATK:           ${params.docker_image_GATK}
- docker_image_DeepSomatic        ${params.docker_image_deepsomatic}
"""

/*
    Nextflow module for converting IntervalList format to BED
    input:
        intervals: path to IntervalList format intervals
        interval_id: interval ID
    params:
        META.workflow_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_GATK: string
*/
process convert_IntervalListToBed_GATK {
    container params.docker_image_GATK

    publishDir path: "${META.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*-contig.bed",
        enabled: params.save_intermediate_files

    ext log_dir: { "${META.log_dir_prefix}/${task.process.split(':')[-1]}/" },
        log_dir_suffix: { "${interval_id}" }

    input:
    val META
    tuple val(interval_id), path(intervals)

    output:
    tuple val(interval_id), path(output_filename), emit: interval_bed

    script:
    output_filename = "${file(intervals).baseName}.bed"
    """
    set -euo pipefail
    gatk IntervalListToBed \
        --INPUT ${intervals} \
        --OUTPUT ${output_filename}
    """
}

/*
    Nextflow module for running DeepSomatic

    input:
        intervals: path to BED format intervals
        interval_id: interval ID
        tumor_bam: path to tumor BAM
        tumor_bam_index: path to tumor BAM index
        normal_bam: path to normal BAM
        normal_bam_index: path to normal BAM index
        reference_fasta: path to reference FASTA
        reference_index: path to reference FASTA index
        reference_dict: path to reference dictionary

    params:
        params.output_dir_base: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_deepsomatic: string
        params.exome: bool.
*/
process call_sSNV_DeepSomatic {
    container params.docker_image_deepsomatic

    tag "${interval_id}"

    publishDir path: "${META.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
               mode: "copy",
               pattern: "*.vcf.gz*",
               enabled: params.save_intermediate_files

    ext log_dir: { "${META.log_dir_prefix}/${task.process.split(':')[-1]}/" },
        log_dir_suffix: { "${interval_id}" }

    input:
    val META
    tuple val(interval_id), path(intervals)
    path(tumor_bam)
    path(tumor_bam_index)
    path(normal_bam)
    path(normal_bam_index)
    path(reference_fasta)
    path(reference_index)
    path(reference_dict)

    output:
    tuple path(vcf_filename), path("${vcf_filename}.tbi"), env(VCF_CALLS), emit: vcf
    tuple path(gvcf_filename), path("${gvcf_filename}.tbi"), env(GVCF_CALLS), emit: gvcf

    script:
    output_filename_base = "${META.output_filename}_unfiltered-${interval_id}"
    vcf_filename = "${output_filename_base}.vcf.gz"
    gvcf_filename = "${output_filename_base}.g.vcf.gz"
    model_type_base = (params.exome) ? "WES" : "WGS"
    sample_mode_extension = (params.sample_mode == 'tumor_only') ? "_TUMOR_ONLY" : ""
    normal_input_args = (params.sample_mode == 'tumor_only') ? "" : "--reads_normal=${normal_bam} --sample_name_normal=\"${params.normal_id}\""
    """
    set -euo pipefail

    mkdir log
    mkdir work

    /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --model_type=${model_type_base}${sample_mode_extension} \
        --ref=${reference_fasta} \
        ${normal_input_args} \
        --reads_tumor=${tumor_bam} \
        --sample_name_tumor="${params.tumor_id}" \
        --output_vcf=${vcf_filename} \
        --output_gvcf=${gvcf_filename} \
        --num_shards=${task.cpus} \
        --logging_dir=log \
        --intermediate_results_dir=work \
        --regions=${intervals} \
        --use_default_pon_filtering=true

    export VCF_CALLS=`zgrep -v ^# ${vcf_filename} | wc -l`
    export GVCF_CALLS=`zgrep -v ^# ${gvcf_filename} | wc -l`
    """
}
