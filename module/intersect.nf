include { compress_file_bzip2; generate_sha512sum } from './common'
include { reorder_samples_BCFtools; intersect_VCFs_BCFtools; plot_VennDiagram_R; concat_VCFs_BCFtools ; convert_VCF_vcf2maf } from './intersect-processes.nf'
include { compress_index_VCF as compress_index_VCF_reordered } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'
include { compress_index_VCF as compress_index_VCF_concat } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'
def sortVcfs(List paths) {
    paths.sort { a, b ->
        def toolA = file(a).getName()
        def toolB = file(b).getName()
        return toolA.compareTo(toolB)
        }
    }
def getToolName(filename) {
    return file(filename).getName().split('-')[0]
    }

workflow intersect {
    take:
    META
    tool_gzvcfs
    tool_indices
    script_dir_ch

    main:
        tool_gzvcfs_ch = tool_gzvcfs
            .flatten()
            .map{ it -> [getToolName(it), it]}
        tool_indices_ch = tool_indices
            .flatten()
        reorder_samples_BCFtools(
            META,
            tool_gzvcfs_ch,
            tool_indices_ch,
            params.tumor_id,
            params.normal_id
            )
        compress_index_VCF_reordered(
            META.combine(reorder_samples_BCFtools.out.gzvcf.map{ it -> ["${getToolName(it)}-SNV", it] })
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
        reordered_indexed_vcfs = compress_index_VCF_reordered.out.index_out
            .map{ it -> [it[0].variant_type, it[1], it[2]] }
        gzvcfs = reordered_indexed_vcfs
            .map{ it -> it[1] }
            .collect()
            .map { sortVcfs(it)  }
        indices = reordered_indexed_vcfs
            .map{ it -> it[2] }
            .collect()
        intersect_VCFs_BCFtools(
            META,
            gzvcfs,
            indices,
            params.intersect_regions,
            params.intersect_regions_index
            )
        plot_VennDiagram_R(
            META,
            script_dir_ch,
            intersect_VCFs_BCFtools.out.isec,
            )
        intersect_vcfs = intersect_VCFs_BCFtools.out.gzvcf
            .flatten()
            .filter{ getToolName(it) != 'DeepSomatic' } // Exclude DeepSomatic for concatenation due to header mis-match
            .collect()
            .map { sortVcfs(it) }
        concat_VCFs_BCFtools(
            META,
            intersect_vcfs,
            intersect_VCFs_BCFtools.out.idx
            )
        convert_VCF_vcf2maf(
            META,
            concat_VCFs_BCFtools.out.vcf,
            params.reference,
            params.normal_id,
            params.tumor_id
            )
        compress_index_VCF_concat(
            META.combine(concat_VCFs_BCFtools.out.vcf.map{ it -> ['SNV', it] })
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
        concat_indexed_vcfs = compress_index_VCF_concat.out.index_out
            .map{ it -> [it[0].variant_type, it[1], it[2]] }
        compress_file_bzip2(META, convert_VCF_vcf2maf.out.maf
            .map{ it -> ['MAF', it]}
            )
        file_for_sha512 = intersect_VCFs_BCFtools.out.gzvcf
            .flatten()
            .map{ it -> ["${getToolName(it)}-vcf", it]}
            .mix(intersect_VCFs_BCFtools.out.idx
                .flatten()
                .map{ it -> ["${getToolName(it)}-idx", it]}
                )
            .mix(concat_indexed_vcfs
                .map{ it -> ["concat-${it[0]}-vcf", it[1]] }
                )
            .mix(concat_indexed_vcfs
                .map{ it -> ["concat-${it[0]}-index", it[2]] }
                )
            .mix(compress_file_bzip2.out.compressed_file
                .map{ it -> ["concat-${it[0]}", it[1]]}
                )
        generate_sha512sum(META, file_for_sha512)
    }
