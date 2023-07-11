
// Estimate quality of raw reads
process run_fastqc {
    conda "${params.conda_envs_path}/fastqc_env"
    tag "$name"
    publishDir "${params.outdir}/reports/fastqc", mode: 'copy'
    cache true
    errorStrategy 'ignore'
 
    input:
        tuple val(name), path(readsForMapping)

    output:
        tuple val(name), path("*.html"), emit: html
        tuple val(name), path("*.zip") , emit: zip

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        fastqc -t 4 -f fastq $readsForMapping
        """
}

// process exon 6 reference
process process_exon6 {
    conda "${params.conda_envs_path}/samtools_and_bwa_env"
    tag "$name"
    publishDir "${params.outdir}/${name.toString()}", mode: 'copy'
    errorStrategy 'ignore'
    cache true
 
    input:
        tuple val(name), path(readsForMapping)
        path(reference)
        path(database)

    output:
        path '*'
        path '*.txt', emit: txt
    
    when:
        task.ext.when == null || task.ext.when
 
    script:
        """
        python $projectDir/bin/AnalyzeAbo_Main.py \\
        --reference="${reference}" \\
        --alleles="${database}" \\
        --output="exon6" \\
        --analysis-type="READS" \\
        --reads="${readsForMapping}" > ${name}_exon6.log.txt 2&>1
        """
}


// process exon 7 reference
process process_exon7 {
    conda "${params.conda_envs_path}/samtools_and_bwa_env"
    tag "$name"
    publishDir "${params.outdir}/${name.toString()}", mode: 'copy'
    errorStrategy 'ignore'
    cache true
 
    input:
        tuple val(name), path(readsForMapping)
        path(reference)
        path(database)

    output:
        path '*'
        path '*.txt', emit: txt
    
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        python $projectDir/bin/AnalyzeAbo_Main.py \\
        --reference="${reference}" \\
        --alleles="${database}" \\
        --output="exon7" \\
        --analysis-type="READS" \\
        --reads="${readsForMapping}" > ${name}_exon7.log.txt 2&>1
        """
}


// compile the tables into single excel
process compile_results{
    tag "Compile tables"
    publishDir "${params.outdir}", mode: 'copy'
 
    input:
        path snp_position_files
        path completed_exon6and7
 
    output:
        path("*.txt")
        path("*.xlsx")
 
    script:
        """
        python $projectDir/bin/Aggregate_ABO_reports.py $snp_position_files
        """
}


// Estimate genome coverage
process bam_stats {
    conda "${params.conda_envs_path}/bedtools_env"
    tag "Coverage and depth: $name"
    publishDir "${params.outdir}/${name}/bam_stats", mode: 'copy'
 
    input:
        tuple val(name), path(mapping_out) 

    output:
        tuple val(name), path("*.coverage"), emit: coverage
        tuple val(name), path("*.average"), emit: average
        tuple val(name), path("*.depth"), emit: depth
        tuple val(name), path("*.consensus"), emit: consensus
        tuple val(name), path("*.stats"), emit: stats
        tuple val(name), path("*.flagstat"), emit: flagstat
        tuple val(name), path("*.tsv"), emit: tsv
    
    script:
        """
        bedtools genomecov -ibam $mapping_out -bg | \\
            awk 'NR>1 {sum[\$1] += \$4; count[\$1] += \$4=="" ? 0 : 1} ; \\
            END {for (i in sum) print i, (count[i] > 0 ? sum[i]/count[i] : "-")}' | \\
            cut -f 1,2 | sort -nrk2 > ${name}.bedtools.coverage

        samtools depth $mapping_out > ${name}.samtools.depth

        awk 'NR>1 {sum[\$1] += \$3; count[\$1] += \$3=="" ? 0 : 1} ; \\
        END {for (i in sum) print i, (count[i] > 0 ? sum[i]/count[i] : "-")}' \\
        ${name}.samtools.depth | \\
        cut -f 1,2 | sort -nrk2 >  ${name}.samtools.depth.average;

        samtools stats $mapping_out | \\
            grep ^SN | \\
            cut -f 2- > ${name}.samtools.consensus

        samtools stats $mapping_out > ${name}.samtools.stats

        samtools flagstat $mapping_out > ${name}.samtools.flagstat

        grep 'error rate' ${name}.samtools.stats | \\
         cut -f 3 | \\
         xargs -I @ echo -e "${name}\t'@'" \\
         > ${name}_consensus_errorrate.tsv
        """
}


// Process MultiQC report
process run_multiqc {
    conda "${params.conda_envs_path}/multiqc_env"
    tag {"MULTIQC Report"}
    publishDir "${params.outdir}/reports/multiqc", mode: 'copy', overwrite: false
    label 'multiqc'
    errorStrategy 'ignore'

    input:
        path multiqc_files

    output:
        path "*multiqc_report.html", emit: report
        path "*_data"              , emit: data
        path "*_plots"             , optional:true, emit: plots

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        multiqc -f $args .
        """
}

// Render Rmarkdown script
process render_rmarkdown {
    tag {"Rmarkdown"}
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
        path(project_dir)
        path(rmd)
        val(db_name)

    output:
        path("*.html")
        path("*.tsv")
        path("*.xlsx")

    script:
        output_f = "${db_name}_raw_read_accuracy.html"
        //path = "$params.outdir/alignments"

        """
        #!/bin/bash -l
        cp -L ${rmd} notebook.Rmd
        grep 'error rate' "${params.outdir}/mapped_reads/bam_stats/*.stats" | cut -f 1,3  > consunsus_error_rate.tsv
        """
       
        """
        #!/usr/bin/env Rscript -e
        ## rm(list=ls())

        rmarkdown::render(
            'render_sequencing_accuracy_report.Rmd', 
            outfile = $output_f, 
            params = list(locus = $db_name), 
            envir = parent.frame()
            )
        """
} 


