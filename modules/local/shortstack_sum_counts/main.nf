process SHORTSTACK_SUM_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d425b12748ce54c44c01a535a1ef5867a6e16cbf62c43151012e893444b1673/data' :
        'community.wave.seqera.io/library/r-base_r-optparse_r-dplyr_r-readr_r-tidyr_bioconductor-rtracklayer:0673f74019be7404' }"

    input:
    tuple val(meta), path(gff), path(counts)

    output:
    path("*.shortstack.counts.vals.txt"), optional: true, emit: counts_txt
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    shortstack_sum_counts.R \\
        --gff $gff \\
        --counts $counts \\
        --count_col $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.shortstack.counts.vals.txt
    """
}
