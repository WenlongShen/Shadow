############################################################
# Dockerfile for shadow
# Shen @ Biological Data Analysis
# Version 0.4
############################################################

# Base image
FROM broadinstitute/gatk:4.6.1.0

# Maintainer
MAINTAINER Wenlong Shen <shenwl1988@gmail.com>

# Install prerequisites
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
        r-base && \
    apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

# Set timezone
RUN ln -sf /usr/share/zoneinfo/Asia/Shanghai /etc/localtime && \
    echo "Asia/Shanghai" > /etc/timezone && \
    dpkg-reconfigure -f noninteractive tzdata
    
# Set up bioconda channels
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# Install packages
RUN conda install conda-libmamba-solver \
                  samtools bedtools deeptools picard fastp \
                  bwa bowtie2 hisat2 stringtie hicexplorer macs3 chip-r -q -y && \
    conda clean -a -q -y
    
