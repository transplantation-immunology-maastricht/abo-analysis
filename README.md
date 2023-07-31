# abo-analysis

Analysis of ABO MinION sequences

ABO sequences were aquired from the NCBI dbRBC database:

<https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home>

# The nextflow workflow

See <https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/> for some literature.

Input files and output directory can be defined in the config files or provided directly in the commandline.

To analyse files with config, run:
  - `nextflow run main.nf -resume ` (user can override inputs and output using `--reads '*.fastq' --outdir 'ABO_results'` on the commandline).

We have also added the ability for the pipeline to automatically set-up a conda or docker based environment with all required tools and libraries. 

Users may also opt for a workload manager such as `-profile slurm,docker|-profile slurm,conda`, is which case, all required modules docker/conda must be installed and loaded. The config slurm parameters must also be defined to ensure tasks are submitted to the correct resource queue/account.

To run without the workload manager but with a specific containerization, use:
  - `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-conda`
  - `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-docker`
