# Base image from ABOTypes Analysis
FROM nvidia/cuda:12.0.0-runtime-ubuntu20.04

# Set the author's name and email as image labels
LABEL image.author.name="Fredrick Mobegi"
LABEL image.author.email="fredrickmaati@gmail.com"

ARG GENERAL_DEPENDENCIES="wget apt-transport-https software-properties-common apt-utils sudo procps cmake"

# Set non-interactive mode to override any user inputs requested by package installations
ENV DEBIAN_FRONTEND noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install --yes --no-install-recommends $GENERAL_DEPENDENCIES && \
    apt-get clean

# Install AzCopy (data movement tool)
RUN mkdir -p /home/azcopy && cd /home/azcopy && \
    wget -O azcopy.tar.gz https://aka.ms/downloadazcopy-v10-linux && \
    tar -xf azcopy.tar.gz --strip=1 
ENV PATH=/home/azcopy:$PATH

# Install PowerShell Core
RUN PLATFORM=$(lsb_release -rs) && \
    wget -q https://packages.microsoft.com/config/ubuntu/${PLATFORM}/packages-microsoft-prod.deb && \
    dpkg -i packages-microsoft-prod.deb && \
    apt-get update && \
    add-apt-repository universe && \
    apt-get install -y powershell

# Install Miniconda
# Specify the path for the installation of conda environments
ARG CONDA_ENVS_PATH="/conda-envs"
# Apply the path to the Conda binaries to the PATH environment variable on the container
ENV PATH="$CONDA_ENVS_PATH/bin:${PATH}"
ARG PATH="$CONDA_ENVS_PATH/bin:${PATH}"
# Download the latest installer for Linux and then run the installation, specifying the path for conda environment
# (-f skips error of folder already existing)
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir -p /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $CONDA_ENVS_PATH \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

# Initialise Conda
RUN conda init

# Finally, set the default path for environments and disable universal activate
# activating base conda breaks cmake gcc+
RUN conda config --append envs_dirs $CONDA_ENVS_PATH
RUN conda config --set auto_activate_base false
RUN conda clean --all --yes

# Install Mamba for faster resolution of Conda packages
RUN conda install -c conda-forge mamba

# Create a Conda environment for Bedtools
RUN mamba create -y -n bedtools_env -c bioconda bedtools

# Create a Conda environment for FastQC
RUN mamba create -y -n fastqc_env -c bioconda fastqc

# Create a Conda environment for MultiQC
RUN mamba create -y -n multiqc_env -c bioconda multiqc

# Create a Conda environment for Samtools and BWA
RUN mamba create -y -n samtools_and_bwa_env -c bioconda samtools bwa minimap2

# Copy the environment.yml file into the container
WORKDIR /home
COPY docker/env.yml env.yml
# Create a new conda environment based on the environment.yml file
RUN mamba create -n abo-analysis-env \
    && mamba env update -n abo-analysis-env -f env.yml \
    && mamba clean --all --yes

RUN echo "source activate abo-analysis-env" >> ~/.bashrc
ENV PATH /conda-envs/abo-analysis-env/bin:$PATH

# Install Python 3 tools
# ARG PYTHON_TOOLS="numpy==1.21.5 Bio==1.5.9 pysam==0.21.0 pandas==1.5.3 scipy matplotlib==3.7.1 biopython openpyxl==3.1.2 xlsxwriter==3.1.2 multiqc==1.17"
ARG PYTHON_TOOLS="numpy==1.26.4 Bio==1.6.2 pysam==0.22.0 pandas==2.2.0 scipy matplotlib==3.8.2 biopython openpyxl==3.1.2 xlsxwriter==3.1.9 multiqc==1.18"

RUN pip3 install --upgrade pip \
    && pip3 install $PYTHON_TOOLS

# Install Nextflow
# Install the latest LTS version of the Adoptium Eclipse Temurin OpenJDK, 
# as recommended for Nextflow (see: https://www.nextflow.io/docs/latest/getstarted.html#requirements)
RUN mkdir -p /etc/apt/keyrings \
    && wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | tee /etc/apt/keyrings/adoptium.asc
RUN echo "deb [signed-by=/etc/apt/keyrings/adoptium.asc] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
RUN apt-get update \
    && apt-get install --yes temurin-17-jdk
# Install Nextflow now that Java is installed
RUN mkdir -p /nextflow_dir && cd /nextflow_dir \
    && wget -qO- https://get.nextflow.io | bash
# Make Nextflow executable for all users, set it to self update, and then add it to PATH
RUN chmod a+rx /nextflow_dir/nextflow
ENV PATH="/nextflow_dir/:${PATH}"
ARG PATH="$/nextflow_dir/:${PATH}"
# Run a first invocation of Nextflow so that it loads all its packages and dependencies
RUN nextflow self-update

ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
# ENV LC_ALL=en_US.UTF-8 <--- This throws an error both when running locally and on Azure
ENV TZ=Australia/Perth

# Define Command or Entry Point
CMD ["/bin/bash"]

# TO BUILD THE CONTAINER RUN:
#    docker build -t pathwestminion.azurecr.io/[**Type Name**]:v[**Type a New Version Number**] -f Dockerfile .  <--Note the period (".") at the end. This is required, and tells Docker to include files in the current folder
# For example:
#    docker build -t pathwestminion.azurecr.io/abotypes1.0:v1 -f Dockerfile .

## -------------------  Additional tests ----------------------------------

# # Use a specific version of micromamba as the base image
# FROM mambaorg/micromamba:0.25.1

# # Set the author's name and email as image labels
# LABEL image.author.name="Fredrick Mobegi"
# LABEL image.author.email="fredrickmaati@gmail.com"

# # Install OpenJDK 11 (Java)
# USER root
# RUN apt-get update && \
#     apt-get install -y --no-install-recommends \
#     openjdk-11-jre wget && \
#     apt-get clean

# # Create a directory for your workflow and set it as the working directory
# WORKDIR /abo-analysis

# # Copy the entire contents of the abo-analysis directory into the container
# COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# # Copy the environment.yml file into the container's /tmp directory
# COPY --chown=$MAMBA_USER:$MAMBA_USER docker/env.yml /tmp/env.yml

# # Create a new conda environment based on the environment.yml file
# RUN micromamba create -n abo-analysis-env -f /tmp/env.yml && \
#     micromamba clean --all --yes

# # Activate the created conda environment
# RUN echo "source activate abo-analysis-env" >> ~/.bashrc
# ENV PATH /opt/conda/envs/abo-analysis-env/bin:$PATH

# # Install Nextflow
# RUN mkdir -p /nextflow_dir && cd /nextflow_dir \
#     && wget -qO- https://get.nextflow.io | bash
# #   Make Nextflow executable for all users, set it to self update, and then add it to PATH
# RUN chmod a+rx /nextflow_dir/nextflow
# ENV PATH="/nextflow_dir/:${PATH}"
# ARG PATH="$/nextflow_dir/:${PATH}"
# #   Run a first invocation of Nextflow so that it loads all its packages and dependencies
# RUN nextflow self-update

# # Install additional system packages
# RUN apt-get update && \
#     apt-get install -y procps g++ && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# # Build
# docker build -t fmobegi/abo-analysis  -f docker/Dockerfile .
