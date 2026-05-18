include { calculate_adjVAF_Python; plot_adjVAF_R } from './plot-vaf-processes'

workflow plot_vaf {
    take:
    META
    identified_gzvcfs
    all_files

    main:
    calculate_adjVAF_Python(
        META,
        identified_gzvcfs,
        all_files
        )

    plot_adjVAF_R(
        META,
        calculate_adjVAF_Python.out.adjusted_vafs
        )
}
