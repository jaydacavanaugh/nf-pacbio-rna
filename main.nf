process get_read_stats {
    tag "fastqc: ${sample}"
    container 'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29'
    
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("*.html"), emit: html
    tuple val(sample), path("*.zip") , emit: zip

    script:
    """
    fastqc -t ${task.cpus} -f bam ${bam}
    """
}

process get_multiqc_stats {
    tag "multiqc"
    container 'community.wave.seqera.io/library/multiqc:1.29--e3ef3b42c5f9f0da'
    publishDir 'results/read_qc', mode: 'copy'

    input:
    path(fqc)
    path(fqz)

    output:
    path("multiqc_report.html")
    path("multiqc_data")

    script:
    """
    multiqc .
    """
}

process cluster_into_isoforms {
    tag "isoseq3: ${sample}"
    container 'community.wave.seqera.io/library/isoseq3:4.0.0--de4f199d3a88eaa8'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.clustered.bam")

    script:
    """
    mkdir -p ${params.scratch}/${sample}
    export TMPDIR=${params.scratch}/${sample} && isoseq3 cluster2 --singletons -j ${task.cpus} "${bam}" ${sample}.clustered.bam
    """
}

process align_sample {
    tag "align: ${sample}"
    container 'community.wave.seqera.io/library/pbmm2:1.17.0--66dffabf917f3cf0'
    publishDir 'results/alignments', mode: 'copy'

    input:
    tuple val(sample), path(bam)
    path(reference)

    output:
    tuple val(sample), path("${sample}.aligned.bam")

    script:
    """
    pbmm2 align --preset ISOSEQ --sort -j ${task.cpus} ${bam} ${reference} ${sample}.aligned.bam
    """
}

process convert_to_sam {
    tag "align: ${sample}"
    container 'aakrosh/bwa-suite:alpine'

    input:
    tuple val(sample), path(bam)
    path(reference)

    output:
    tuple val(sample), path("${sample}.sam")

    script:
    """
    samtools calmd -@ ${task.cpus} --output-fmt SAM ${bam} ${reference} > ${sample}.sam
    """
}

process label_reads_using_talon {
    tag "label: ${sample}"
    container 'aakrosh/talon:v6'
    publishDir 'results/talon', mode: 'copy'

    input:
    tuple val(sample), path(sam)
    path(reference)
    path(faidx)

    output:
    tuple val(sample), path("${sample}_labeled.sam"), path("${sample}_read_labels.tsv")    

    script:
    """
    talon_label_reads --f ${sam} --g ${reference} --t ${task.cpus} \
        --tmpDir ${params.scratch}/${sample} --deleteTmp \
        --o ${sample}
    """    
}

process initialize_run_talon {
    tag 'TALON db'
    container 'aakrosh/talon:v6'
    publishDir 'results/talon', mode: 'copy'

    input:
    path(gtf)
    val(sample_info_list)

    output:
    path("talon.db"), emit: talon_db
    tuple path("talon_swan_QC.log"), path("talon_swan_talon_read_annot.tsv"), emit: annot

    script:
    def lines = sample_info_list.collect { tuple ->
        def sample_id = tuple[0]
        def sam_file = file(tuple[1]).toAbsolutePath()
        return "${sample_id},${sample_id},${params.platform},${sam_file}"
    }.join('\n')

    """
    talon_initialize_database --f ${gtf} \
    --g ${params.genome} --a gencode --l 0 --idprefix TALON --5p 500 --3p 300 \
    --o talon

    echo "${lines}" > config.csv
 
    talon --f config.csv --db talon.db -t ${task.cpus} \
        --build ${params.genome} --cov 0.9 --identity 0.8 --o talon_swan
    """
}

process filter_transcripts {
    tag 'TALON filter'
    container 'aakrosh/talon:v6'

    input:
    path(db)

    output:
    path("pass_list.csv")

    script:
    """
    talon_filter_transcripts --db ${db} -a gencode --maxFracA 0.5 \
        --minCount 2 --minDatasets 2 --o pass_list.csv
    """
}

process create_gtf {
    tag 'TALON gtf'
    container 'aakrosh/talon:v6'
    publishDir 'results/talon', mode: 'copy'
    
    input:
    path(pass_file)
    path(db)

    output:
    path("all_talon_observedOnly.gtf")

    script:
    """
    talon_create_GTF --db ${db} -b ${params.genome} -a gencode \
      --whitelist ${pass_file} --observed --o all
    """
}

process get_abundance {
    tag 'TALON abundance'
    container 'aakrosh/talon:v6'
    publishDir 'results/talon', mode: 'copy'

    input:
    path(pass_file)
    path(db)
  
    output:
    path("all_talon_abundance_filtered.tsv")

    script:
    """
    talon_abundance --db ${db} -a gencode -b ${params.genome} \
        --whitelist ${pass_file} --o all
    """
}

process get_isoform_stats {
    tag 'squanti'
    container 'anaconesalab/sqanti3:v5.4'
    containerOptions { "--cleanenv" }
    publishDir 'results/isoform_qc', mode: 'copy'

    input:
    path(talon_gtf)
    path(gtf)
    path(reference)

    output:
    path("isoforms_*")

    script:
    """
    conda run --no-capture-output -n sqanti3 sqanti3_qc.py -t ${task.cpus} --force_id_ignore -o isoforms --report html ${talon_gtf} ${gtf} ${reference}
    """
}

workflow {
    if (!params.apptainer_cache_dir) {
        error "Missing required parameter: --apptainer_cache_dir"
    }
    if (!params.queue) {
       error "Missing required parameter: --queue" 
    }
    if (!params.clusterOptions) {
       error "Missing required parameter: --clusterOptions"
    }
    if (!params.flnc_bams) {
       error "Missing required parameter: --flnc_bams"
    }
    if (!params.gtf) {
        error "Missing required parameter: --gtf"
    }
    if (!params.reference) {
        error "Missing required parameter: --reference"
    }
    if (!params.reference_faidx) {
        error "Missing required parameter: --reference_faidx"
    }
    if (!params.scratch) {
        error "Missing required parameter: --scratch"
    }
    if (!params.genome) {
        error "Missing required parameter: --genome"
    }
    if (!params.platform) {
        error "Missing required parameter: --platform"
    }
    
    Channel
        .from(params.flnc_bams)
        .map { [ it.id, it.path ] }
        .set { input_channel }

    // run QC on the FLNC BAMs
    get_read_stats(input_channel)
        .set { qc_channel }

    qc_channel
        .html.map { it[1] }
        .collect()
        .set { all_html_channel }

    qc_channel
        .zip.map { it[1] }
        .collect()
        .set { all_zip_channel }


    get_multiqc_stats(all_html_channel, all_zip_channel)

    // Cluster FLNC BAM into Isoforms
    cluster_into_isoforms(input_channel)
        .set { clustered_channel }

    // align using minimap2
    align_sample(clustered_channel, file(params.reference))
        .set { aligned_channel }

    // convert the BAM to SAM
    convert_to_sam(aligned_channel, file(params.reference))
        .set { sam_channel }

    // Label reads with TALON for detection of internal priming
    label_reads_using_talon(sam_channel, file(params.reference), file(params.reference_faidx))
        .set { label_channel }

    // initialize and run TALON 
    label_channel
        .collect(flat: false)
        .set { collected_label_channel }

    initialize_run_talon(file(params.gtf), collected_label_channel)
        .set { talon_channel }
    
    // Filter annotated transcripts for internal priming and reproducibility
    filter_transcripts(talon_channel.talon_db)
        .set { talon_filtered_channel }

    // Create GTF of transcripts for observed transcripts after filtering
    create_gtf(talon_filtered_channel, talon_channel.talon_db)
        .set { gtf_channel }

    // Obtain a filtered abundance file of each transcript
    get_abundance(talon_filtered_channel, talon_channel.talon_db)
        .set { abundance_channel }    

    // get some isoform stats
    get_isoform_stats(gtf_channel, file(params.gtf), file(params.reference))
}
