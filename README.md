# ABO blood typing using Oxford Nanopore MinION sequencing.

ABO sequences were aquired from the NCBI dbRBC database:

[https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home](https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home)

See [https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/](https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/) for some literature.

# Required tools

The pipeline makes use of the following dependencies:

- bioconda::fastqc=0.11.9
- bioconda::bwa=0.7.17
- bioconda::samtools=1.16.1
- bioconda::multiqc=1.14
- conda-forge::biopython=1.81
- python=3.10
- pip:
  - numpy==1.21.5
  - Bio==1.5.9
  - openpyxl==3.1.2
  - pandas==1.5.3
  - pysam==0.21.0
  - matplotlib==3.7.1
  - xlsxwriter==3.1.2

# Testing without `nextflow`

The pipeline can be tested of single input file by cloning this repo and installing all dependncies above, then running the following commands:

```
python bin/AnalyzeAbo_Main.py  \
 --reference="assets/A1_01_01_1_reference_Exon6.fasta" \
 --alleles="assets/ABO_Database.fasta" \
 --output="SampleName/exon6" \
 --analysis-type="READS" \
 --reads="SampleName.fastq" \
 
 python bin/AnalyzeAbo_Main.py  \
 --reference="assets/reads_bc51/A1_01_01_1_reference_Exon7.fasta" \
 --alleles="assets/input/ABO_Database.fasta" \
 --output="SampleName/exon7" \
 --analysis-type="READS" \
 --reads="SampleName.fastq"
```

Looping through a couple of samples with the above command will generate the following outputs:

Data structure

```
OutputDirectoryName/
├── Sample1
│   ├── exon6
│   │   └── alignment
│   └── exon7
│       └── alignment
├── Sample2
│   ├── exon6
│   │   └── alignment
│   └── exon7
│       └── alignment
├── Sample3
│   ├── exon6
│   │   └── alignment
│   └── exon7
│       └── alignment

```

With individual files named as follows:

```
OutputDirectoryName/
├── Sample1
│   ├── exon6
│   │   ├── ABOPhenotype.txt
│   │   ├── ABOReadPolymorphisms.txt
│   │   ├── alignment
│   │   │   ├── alignment.bam
│   │   │   ├── alignment.bam.bai
│   │   │   ├── AlignmentReference.fasta
│   │   │   ├── AlignmentReference.fasta.amb
│   │   │   ├── AlignmentReference.fasta.ann
│   │   │   ├── AlignmentReference.fasta.bwt
│   │   │   ├── AlignmentReference.fasta.pac
│   │   │   └── AlignmentReference.fasta.sa
│   │   └── ReadAlignmentSpreadsheet.csv
│   ├── exon7
│   │   ├── ABOPhenotype.txt
│   │   ├── ABOReadPolymorphisms.txt
│   │   ├── alignment
│   │   │   ├── alignment.bam
│   │   │   ├── alignment.bam.bai
│   │   │   ├── AlignmentReference.fasta
│   │   │   ├── AlignmentReference.fasta.amb
│   │   │   ├── AlignmentReference.fasta.ann
│   │   │   ├── AlignmentReference.fasta.bwt
│   │   │   ├── AlignmentReference.fasta.pac
│   │   │   └── AlignmentReference.fasta.sa
│   │   └── ReadAlignmentSpreadsheet.csv
│   ├── Sample1_exon6.log.txt
│   └── Sample1_exon7.log.txt

```

The `ABOPhenotype.txt` files from each sampe can then be collated using:

`python bin/Aggregate_ABO_reports.py OutputDirectoryName`

# The `nextflow` workflow

The steps above are simplified in a `NextFlow; https://www.nextflow.io/` pipeline that does all the above steps and streamlines installation of requisite software and tools with a single command.

Besides reproducability, nextflow offeres several advatages over conventional `for loops`, including scallability, portability, and debugging/resumption of failed tasks.

Input files and output directory can be defined in the config files or provided directly in the commandline.

To analyse files with config, run:

- `nextflow run main.nf -resume ` (user can override inputs and output using `--reads '*.fastq' --outdir 'ABO_results'` on the commandline).

We have also added the ability for the pipeline to automatically set-up a conda or docker based environment with all required tools and libraries.

Users may also opt for a workload manager such as `-profile slurm,docker|-profile slurm,conda`, is which case, all required modules docker/conda must be installed and loaded. The config slurm parameters must also be defined to ensure tasks are submitted to the correct resource queue/account.

To run without the workload manager but with a specific containerization, use:

- `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-conda` or `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -profile conda`
- `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-docker` or `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -profile docker`

# Results from the `Nextflow` pipeline will look something like this:

```
230128R_ABO_results/
├── ABO_result.txt
├── ABO_result.xlsx
├── execution_report.html
├── execution_timeline.html
├── execution_trace.txt
├── SampleName
│   ├── exon6
│   │   ├── ABOPhenotype.txt
│   │   ├── ABOReadPolymorphisms.txt
│   │   ├── alignment
│   │   │   ├── alignment.bam
│   │   │   ├── alignment.bam.bai
│   │   │   ├── AlignmentReference.fasta
│   │   │   ├── AlignmentReference.fasta.amb
│   │   │   ├── AlignmentReference.fasta.ann
│   │   │   ├── AlignmentReference.fasta.bwt
│   │   │   ├── AlignmentReference.fasta.pac
│   │   │   └── AlignmentReference.fasta.sa
│   │   └── ReadAlignmentSpreadsheet.csv
│   ├── exon7
│   │   ├── ABOPhenotype.txt
│   │   ├── ABOReadPolymorphisms.txt
│   │   ├── alignment
│   │   │   ├── alignment.bam
│   │   │   ├── alignment.bam.bai
│   │   │   ├── AlignmentReference.fasta
│   │   │   ├── AlignmentReference.fasta.amb
│   │   │   ├── AlignmentReference.fasta.ann
│   │   │   ├── AlignmentReference.fasta.bwt
│   │   │   ├── AlignmentReference.fasta.pac
│   │   │   └── AlignmentReference.fasta.sa
│   │   └── ReadAlignmentSpreadsheet.csv
│   ├── SampleName_exon6.log.txt
│   └── SampleName_exon7.log.txt
├── software_versions.txt
└── workflow.oncomplete.txt
```

Feel free to reach out if you need any support getting this tool running or with suggestions for improvement.
