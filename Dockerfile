# nf-rt-rnaseq task container — bio tools + the FDZ004.2 pipeline scripts. Build for linux/amd64 (AWS Batch).
FROM mambaorg/micromamba:1.5.8
USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps curl unzip ca-certificates && \
    rm -rf /var/lib/apt/lists/*
# AWS CLI v2 — used inside the task for the Secrets Manager fetch (S3 staging is handled by Nextflow)
RUN curl -sSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip && \
    unzip -q /tmp/awscliv2.zip -d /tmp && /tmp/aws/install && rm -rf /tmp/aws /tmp/awscliv2.zip
# minimal bioinformatics core (matches the pipeline spec; benchmark aligners dropped)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.12 'sra-tools>=3.0' bowtie2 minimap2 samtools bedtools pigz \
      python-duckdb pysam numpy requests && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
# vendored pipeline code (see README: copy the needed scripts into ./scripts before building)
COPY scripts/ /pipeline/scripts/
WORKDIR /work
