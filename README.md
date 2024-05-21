# ABO blood typing using Oxford Nanopore MinION sequencing

ABO sequences were aquired from the NCBI dbRBC database:

[https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home](https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home)

See [https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/](https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/) for some literature.

# Required tools

The pipeline makes use of the following core dependencies:

- bioconda::fastqc=0.12.1
- bioconda::bwa=0.7.17
- conda-forge::ncurses
- bioconda::samtools=1.19.2
- bioconda::minimap2=2.26
- conda-forge::biopython=1.83
- python=3.10
- pip
- pip:
  - numpy>=1.26.0
  - Bio>=1.6.0
  - biopython>=1.8o
  - openpyxl>=3.1.0
  - pandas>=2.2.0
  - pysam>=0.22.0
  - matplotlib>=3.8.0
  - XlsxWriter>=3.2.0
  - multiqc>=1.18

A complete list of dependencies is found in the assets folder `assets/conda.yml`.

# Testing without `nextflow`

The pipeline can be tested of single input file by cloning this repo and installing all dependncies above, then running the following commands:

```python
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

```bash
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

```bash
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

- `nextflow run main.nf -resume` (user can override inputs and output using `--reads '*.fastq' --outdir 'ABO_results'` on the commandline).

We have also added the ability for the pipeline to automatically set-up a conda or docker based environment with all required tools and libraries.

Users may also opt for a workload manager such as `-profile slurm,docker|-profile slurm,conda`, is which case, all required modules docker/conda must be installed and loaded. The config slurm parameters must also be defined to ensure tasks are submitted to the correct resource queue/account.

For `conda` environment, it is advisable to prepare the working computer using mamba for easy resolution of environments.
Follow these steps to achieve better results.

```bash
mamba create -y -n abo-analysis-env
conda activate abo-analysis-env
mamba env update --file abo-analysis/assets/conda.yml --prune
conda deactivate

# If conda cativate fails, run:
source {path_to_anaconda}/anaconda3/etc/profile.d/conda.sh
```

To run without the workload manager but with a specific containerization, use:

- `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-conda abo-analysis-env` or
  `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -profile conda`
- `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-docker fmobegi/abo-analysis` or
  `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -profile docker`

# Renaming samples (OPTIONAL)
**NB: This step is on by default. If not needed, users must explicitly disable it in the config file by setting** `--skip_renaming true`

The code by default renames samples using a tab file with `sequencingID` and `sampleName` (see `nextflow.config` file under `$params.renaming_file`).
In our use case, this file is exported directly from LIS Soft into an excel file whose "Acc#" and "Patient Name" columns (A,C) are used for renaming purposes. 
The script `rename_samples.py` can be modified according to fit each use case. 
This option is controlled by the parameter `$params.skip_renaming` and can be overridden via the commandline using option `--skip_renaming true` to skip the process.
Please feel free to reach our if you need help customising this part for your use-case

# Results from the `Nextflow` pipeline will look something like this

```bash
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

Feel free to raise an issue or reach out if you need any support getting this tool running, or with suggestions for improvement.
