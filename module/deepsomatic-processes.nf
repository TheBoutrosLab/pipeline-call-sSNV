log.info """\
====================================
          M U T E C T 2
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
        params.workflow_output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
*/
process convert_IntervalListToBed_GATK {
    container params.docker_image_GATK

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*-contig.bed",
        enabled: params.save_intermediate_files

    ext log_dir: { "DeepSomatic-${params.deepsomatic_version}/${task.process.split(':')[-1]}" }
    ext log_dir_suffix: { "${interval_id}" }

    input:
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
    Nextflow module for running DeepVariant

    input:
        sample_id: sample ID
        bam: path to sample BAM
        bam_index: path to BAM index
        intervals: path to BED format intervals
        interval_id: interval ID

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_deepvariant: string
*/
process call_gSNP_DeepVariant {
    container params.docker_image_deepsomatic

    tag "${sample_id}-${interval_id}"

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':')[-1]}",
               mode: "copy",
               pattern: "*.vcf.gz*",
               enabled: params.save_intermediate_files

    ext log_dir: { "DeepSomatic-${params.deepsomatic_version}/${task.process.split(':')[-1]}" }
    ext log_dir_suffix: { "${interval_id}" }

    input:
    path(normal_bam)
    path(normal_bam_index)
    path(tumor_bam)
    path(tumor_bam_index)
    tuple path(intervals), val(interval_id)
    path(reference_fasta)
    path(reference_index)
    path(reference_dict)

    output:
    tuple path(vcf_filename), path("${vcf_filename}.tbi"), env(VCF_CALLS), emit: vcf
    tuple path(gvcf_filename), path("${gvcf_filename}.tbi"), env(GVCF_CALLS), emit: gvcf

    script:
    output_filename_base = "${params.output_filename}_unfiltered-${interval_id}"
    vcf_filename = "${output_filename_base}.vcf.gz"
    gvcf_filename = "${output_filename_base}.g.vcf.gz"
    model_type = (params.exome) ? "WES" : "WGS"
    """
    set -euo pipefail

    mkdir log
    mkdir work

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=${model_type} \
        --ref=${reference_fasta} \
        --reads_normal=${normal_bam} \
        --reads_tumor=${tumor_bam} \
        --sample_name_normal="${params.normal_id}" \
        --sample_name_tumor="${params.tumor_id}" \
        --output_vcf=${vcf_filename} \
        --output_gvcf=${gvcf_filename} \
        --num_shards=${task.cpus} \
        --logging_dir=log \
        --intermediate_results_dir=work \
        --regions=${intervals} \

    export VCF_CALLS=`zgrep -v ^# ${vcf_filename} | wc -l`
    export GVCF_CALLS=`zgrep -v ^# ${gvcf_filename} | wc -l`
    """
}
