# ABO blood typing using Oxford Nanopore MinION sequencing

ABO sequences were aquired from the NCBI dbRBC database:

[https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home](https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home).

The data is currently archived at https://ftp.ncbi.nlm.nih.gov/pub/mhc/mhc/Final%20Archive/

See [https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/](https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/) for some literature.

## Required tools

The pipeline makes use of the following core dependencies:

```yaml
- bioconda::fastqc=0.12.1
- bioconda::bwa=0.7.17
- conda-forge::ncurses=6.4.20240210
- bioconda::samtools=1.19.2
- bioconda::minimap2=2.27
- conda-forge::biopython=1.83
- python=3.10
- pip
- pip:
  - numpy>=1.26.0
  - Bio>=1.6.2
  - biopython>=1.83
  - openpyxl>=3.1.9
  - pandas>=2.2.1
  - pysam>=0.22.0
  - matplotlib>=3.8.3
  - XlsxWriter>=3.2.0
  - multiqc>=1.21
```

A complete list of dependencies is found in the assets folder `assets/conda.yml`.

## Testing without `nextflow`

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

```yaml
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

```yaml
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

## The `nextflow` workflow

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

# If conda activate fails, run:
source {path_to_anaconda}/anaconda3/etc/profile.d/conda.sh
```

## Some housekeeping for general usage

This pipeline was optimized to run our in-house samples whose naming convention may be different from other users.
The in-house sequencing naming convention starts with the alphabetical initials `IMM` (result) or `INGS` (+ve control) or `NGS` (-ve control), followed by a dash `-` and a 2-digit year code `##`, followed by another dash `-` and a numeric sample id `####+` for the sequencing run, followed by an `underscore (_)` and the barcode number `barcode##`. Spliting of these sequencing IDs into `samplenale` and `barcode##`, required to aggregate the results, is performed by the `process_files` function in the `Aggregate_ABO_reports.py` script.
This function needs to be modified to suit individual needs based on the sequencing IDs in your sequencing core.

This `pipeline` can also process `short-reads` from both `Illumina` and `IonTorrent` platforms as long as the paired reads are interleaved into a single file per sample and the naming convention follows the `IMM-##-####+_barcode##` convention.

### `process_files` function

```python
def process_files(self):
    for filename in os.listdir(self.input_dir):
        if os.path.isdir(os.path.join(self.input_dir, filename)):
            # match = re.match(r"^IMM-[0-9]+-[0-9]+_barcode\d+", filename)
            match = re.match(r"^(IMM|INGS|NGS)(-[0-9]+-[0-9]+)?_barcode\d+", filename)
            if match:
                print("Processing file: " + filename)
                # Extract barcode and sample_name from the filename
                sample_name, barcode = filename.split("_")
                self.process_file(filename)
                print("Done adding Sample %s with barcode %s to merged data frame\n" % (sample_name, barcode))

```

To run without the workload manager but with a specific containerization setup, use:

### For conda

```bash
nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-conda abo-analysis-env
```

or

```bash
nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -profile conda
```

and

### For docker

```bash
nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-docker fmobegi/abo-analysis
```

or

```bash
nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -profile docker
```

## Results

The `Nextflow` pipeline output directory structure will look something like this:

```yaml
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

## (Un)happy with the pipeline?

Feel free to leave a comment in the discussions section `https://github.com/fmobegi/abo-analysis/discussions/2` to let us know your experience using this tool, or raise an issue (<https://github.com/transplantation-immunology-maastricht/abo-analysis/issues>) or reach out if you need any support getting this tool running, or with suggestions for improvement.
