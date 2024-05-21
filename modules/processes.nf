
// Estimate quality of raw reads
process run_fastqc {
    tag "$name"
    label 'process_medium'

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
        fastqc -t $task.cpus -f fastq $readsForMapping
        """
}


// process exon 6 reference
process process_exon6 {
    tag "$name"
    label 'process_high'

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
        python3 $projectDir/bin/AnalyzeAbo_Main.py \\
            --reference="${reference}" \\
            --alleles="${database}" \\
            --output="exon6" \\
            --analysis-type="READS" \\
            --reads="${readsForMapping}" 2>&1 | tee  ${name}_exon6.log.txt
        """
}


// process exon 7 reference
process process_exon7 {
    tag "$name"
    label 'process_high'

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
        python3 $projectDir/bin/AnalyzeAbo_Main.py \\
            --reference="${reference}" \\
            --alleles="${database}" \\
            --output="exon7" \\
            --analysis-type="READS" \\
            --reads="${readsForMapping}" 2>&1 | tee ${name}_exon7.log.txt
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
        path("*.csv")
        path("*.log")
        path("final_export.csv"),  emit: final_export
 
    script:
        """
        python3 $projectDir/bin/Aggregate_ABO_reports.py $snp_position_files 2>&1 | tee  ABO_results.log
        """
}


// Reformat samples 
process rename_samples {
    tag "Handle Final Export"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(deobfuscation)
    path(final_export_file)
    

    output:
    path "MatchPointExport.txt"
    path "MatchPointExport_with_sequencingAcc.txt"

    script:
    """
    python3 $projectDir/bin/rename_samples.py ${final_export_file} ${deobfuscation}
    """
}


// Process MultiQC report
process run_multiqc {
    tag {"MULTIQC Report"}
    publishDir "${params.outdir}/reports/multiqc", mode: 'copy', overwrite: false
    label 'process_single'

    errorStrategy 'ignore'

    input:
        // path outdir
        path multiqc_files
        path logo

    output:
        path "*_multiqc.html",      emit: report
        path "*_multiqc_data",      emit: data
        path "*_multiqc_plots",     optional:true, emit: plots

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        python3 $projectDir/bin/generate_samples_mapping.py ${params.outdir} 2>&1 | tee generate_samples_mapping.log 
        
        multiqc -ip -f --ignore "alignment.*"  --filename ABO_typing_multiqc $args .
        rm *.samtools.*
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
