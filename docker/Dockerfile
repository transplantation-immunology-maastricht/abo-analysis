# Use a specific version of micromamba as the base image
FROM mambaorg/micromamba:0.25.1

# Set the author's name and email as image labels
LABEL image.author.name="Fredrick Mobegi"
LABEL image.author.email="fredrickmaati@gmail.com"

# Copy the environment.yml file into the container's /tmp directory
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

# Create a new conda environment based on the environment.yml file
USER root
RUN micromamba create -n abo-analysis-env -f /tmp/env.yml && \
    micromamba clean --all --yes 

# Activate the created conda environment and install additional packages from environment.yml
RUN echo "source activate abo-analysis-env" >> ~/.bashrc
ENV PATH /opt/conda/envs/abo-analysis-env/bin:$PATH

# Install additional system packages
RUN apt-get update && \
    apt-get install -y procps g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Clean up the package manager caches
USER $MAMBA_USER

## Build
#  docker build -t fmobegi/abo-analysis .
## Push
#