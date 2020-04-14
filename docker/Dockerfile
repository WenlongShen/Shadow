############################################################
# Dockerfile for Shadow
# Shen @ data analysis on genome-wide studies
# Version 0.2
############################################################

# Base image
FROM broadinstitute/gatk:4.1.5.0

# Maintainer
MAINTAINER Wenlong Shen <shenwl1988@gmail.com>

# Install prerequisites
#RUN apt-get update && apt-get upgrade -y && apt-get install -y \
#        r-base && \
#    apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

# Set timezone
RUN ln -sf /usr/share/zoneinfo/Asia/Shanghai /etc/localtime && \
    echo "Asia/Shanghai" > /etc/timezone && \
    dpkg-reconfigure -f noninteractive tzdata
    
# Set up bioconda channels
RUN conda config --add channels bioconda

# Install packages
RUN conda install fastp bwa -q -y && \
    conda clean -a -q -y

# Install Illumina manta & strelka
RUN cd /home/ && \
    wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 && \
    tar xvjf strelka-2.9.10.centos6_x86_64.tar.bz2 && \
    wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 && \
    tar xvjf manta-1.6.0.centos6_x86_64.tar.bz2 && \
    rm *.bz2

ENV PATH="${PATH}:/home/strelka-2.9.10.centos6_x86_64/bin:/home/manta-1.6.0.centos6_x86_64/bin"
