process DESEQ2_QC {
    label "process_medium"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d425b12748ce54c44c01a535a1ef5867a6e16cbf62c43151012e893444b1673/data' :
        'community.wave.seqera.io/library/r-base_r-optparse_r-ggplot2_r-rcolorbrewer_pruned:9e75394d0bc21987' }"

    input:
    tuple val(meta), path(counts)

    output:
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*rle.vals.txt"        , optional:true, emit: rle_txt
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('bioconductor-deseq2'), eval("Rscript -e \"library(DESeq2); cat(as.character(packageVersion('DESeq2')))\""), emit: versions_deseq2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    def label_upper = args2.toUpperCase()
    prefix = task.ext.prefix ?: "deseq2"
    """
    deseq2_qc.R \\
        --count_file $counts \\
        --outdir ./ \\
        --cores $task.cpus \\
        --outprefix $prefix \\
        $args
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    prefix = task.ext.prefix ?: "deseq2"
    """
    touch ${prefix}.dds.RData
    touch ${prefix}.pca.vals.txt
    touch ${prefix}.rle.vals.txt
    touch ${prefix}.plots.pdf
    touch ${prefix}.sample.dists.txt
    touch R_sessionInfo.log

    mkdir size_factors
    touch size_factors/${prefix}.size_factors.RData
    for i in `head $counts -n 1 | cut -f3-`;
    do
        touch size_factors/\${i}.size_factors.RData
    done
    """
}
