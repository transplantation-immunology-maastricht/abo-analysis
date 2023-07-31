// Find the versions of each required tool.
// Will fail if any of the required software are missing.

workflow check_env {
    main:
    fastqc = get_fastqc_version()
    bwa = get_bwa_version()
    samtools = get_samtools_version()
    multiqc = get_multiqc_version()
    
	emit:
    fastqc
    bwa
    samtools
    multiqc
}

process get_fastqc_version {

    label 'fastqc'

    output:
    env VERSION

    script:
    """
    if ! which fastqc > /dev/null
    then
        echo -e "Could not find the program 'fastqc' in your environment path.\n" 1>&2

        echo "Please install fastqc before you can continue." 1>&2

        exit 127
    fi

    VERSION="\$(fastqc --version)"
	
    """
}

process get_bwa_version {

    label 'BWA'

    output:
    env VERSION

    script:
    """
    #!/bin/bash

    if ! which bwa > /dev/null; then
        echo "Could not find the program 'bwa' in your environment path." >&2
        echo "Please install BWA before you can continue." >&2
        exit 127
    fi

    VERSION="\$(bwa 2>&1 | grep -oP 'Version: \\K\\S+')"
    """
}

process get_samtools_version {

    label 'samtools'
	label 'htslib'

    output:
    env VERSION

    script:
    """
    if ! which samtools > /dev/null
    then
        echo -e "Could not find the program 'samtools' in your environment path.\n" 1>&2

        echo "Please install samtools before you can continue." 1>&2

        exit 127
    fi

    VERSION="\$(samtools --version| grep -E "(^samtool|Using htslib)"  | sed -e "s/Using/: using/g" )"
    """
}

process get_multiqc_version {
    label 'multiqc'

    output:
    env VERSION

    script:
    """
    if ! which multiqc > /dev/null
    then
        echo -e "Could not find the program 'multiqc' in your environment path.\n" 1>&2

        echo "Please install multiqc before you can continue." 1>&2

        exit 127
    fi

    VERSION="\$( multiqc --version )"
    """
}

process publish_software {
  tag 'Software_versions'

  publishDir params.outdir, mode: 'copy'

  output:
    path 'software_versions.txt', emit: txt

  script:
  """
    echo "\n---------------------------------------" > tmp.txt
    echo "Software used for this analysis
    \n---------------------------------------" >> tmp.txt
    echo "
    FASTQC: "\$(fastqc --version | sed -e "s/FastQC v//g")"
    BWA: "\$(bwa 2>&1 | grep -oP 'Version: \\K\\S+')"
    SAMTOOLS:"\$(samtools --version| grep -E "(^samtool|Using htslib)"  | \\
        sed -e "s/samtools//g" |\\
        sed -e "s/Using/: using/g" )""  >> tmp.txt
    echo "
    MULTIQC: "\$( multiqc --version | sed -e "s/multiqc, version //g")" 
    " >> tmp.txt
    echo "\n---------------------------------------" >> tmp.txt
    sed -i 's/^ *//g'  tmp.txt
    sed -i '/^\$/d'  tmp.txt
    cp tmp.txt software_versions.txt
    rm tmp.txt
  """
}

