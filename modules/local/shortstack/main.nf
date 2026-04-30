process SHORTSTACK {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/shortstack:4.1.0--hdfd78af_0'
        : 'quay.io/biocontainers/shortstack:4.1.0--hdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    path  fasta

    output:
    tuple val(meta), path("*merged_condensed.bam")     , optional: true, emit: bam
    tuple val(meta), path("*merged_condensed.bam.csi") , optional: true, emit: csi
    tuple val(meta), path("*merged_alignments.bam")    , optional: true, emit: alignment_bam
    tuple val(meta), path("*merged_alignments.bam.csi"), optional: true, emit: alignment_csi
    tuple val(meta), path("*.counts.txt")              , optional: true, emit: counts
    tuple val(meta), path("*.gff3")                    , emit: gff3
    tuple val(meta), path("*.results.txt")             , emit: results
    tuple val(meta), path("*.txt")                     , emit: txt
    tuple val("${task.process}"), val('shortstack'), eval('ShortStack --version | sed "s/ShortStack //"'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def genome_fasta = fasta.toString().endsWith(".fasta") || fasta.toString().endsWith(".fa") ? fasta : "${fasta.baseName}.fasta"
    def is_bam       = reads.toString().endsWith(".bam")
    def input_arg    = is_bam ? "--bamfile $reads" : "--readfile $reads"
    output_dir_name  = "${prefix}"
    """
    if [ ! -f ${genome_fasta} ]; then
        ln -sf $fasta ${genome_fasta}
    fi

    ShortStack \\
        --genomefile ${genome_fasta} \\
        ${input_arg} \\
        --threads $task.cpus \\
        --outdir $prefix \\
        $args

    if [ -f ${prefix}/Results.gff3 ]; then
        mv ${prefix}/Results.gff3 ${prefix}/${prefix}.gff3
    fi
    if [ -f ${prefix}/Results.txt ]; then
        mv ${prefix}/Results.txt ${prefix}/${prefix}.results.txt
    fi
    if [ -f ${prefix}/Counts.txt ]; then
        mv ${prefix}/Counts.txt ${prefix}/${prefix}.counts.txt
    fi

    mv ${prefix}/* .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_dir_name = "${prefix}"
    """
    mkdir -p ${prefix}
    touch ${prefix}.bam
    touch ${prefix}.bam.csi
    touch ${prefix}.gff3
    touch ${prefix}.results.txt
    touch ${prefix}.counts.txt
    """
}