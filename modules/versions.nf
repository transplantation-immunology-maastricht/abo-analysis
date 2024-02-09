// Find the versions of each required tool.
// Will fail if any of the required software are missing.

workflow check_env {
    main:
        fastqc = get_fastqc_version()
        bwa = get_bwa_version()
        minimap2 = get_minimap2_version()
        samtools = get_samtools_version()
        multiqc = get_multiqc_version()
    
	emit:
        fastqc
        bwa
        minimap2
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

process get_minimap2_version {
    label 'minimap2'

    output:
    env VERSION

    script:
    """
    if ! which minimap2 > /dev/null
    then
        echo -e "Could not find the program 'minimap2' in your environment path.\n" 1>&2

        echo "Please install minimap2 before you can continue." 1>&2

        exit 127
    fi

    VERSION="\$(minimap2 --version)"
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
    path 'software_mqc_versions.yml', emit: txt

  script:
  """
  echo "\n" > tmp.txt
  echo "
  Raw reads QC:
      fastqc: \'"\$(fastqc --version | sed -e "s/FastQC v//g")"\'
  Mapping and BAM processing:
      bwa: \'"\$(bwa 2>&1 | grep -oP 'Version: \\K\\S+')"\'
      minimap2: \'"\$(minimap2 --version)"\'
      samtools:
          samtools: \'"\$(samtools --version| grep -E "(^samtool)" | sed -e "s/samtools //g")"\'
          htslib: \'"\$(samtools --version| grep -E "(^Using htslib)" | sed -e "s/Using htslib //g")"\'
      "  >> tmp.txt
  echo "
  Reporting:
      multiqc: \'"\$( multiqc --version | sed -e "s/multiqc, version //g")"\'
  " >> tmp.txt
  sed -i '/^\$/d'  tmp.txt
  cp tmp.txt software_mqc_versions.yml
  rm tmp.txt
  """
}
