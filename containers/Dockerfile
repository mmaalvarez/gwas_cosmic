# Use Ubuntu Jammy as the base image
FROM ubuntu:jammy

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C

# Install initial packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    make \
    gcc \
    g++ \
    gfortran \
    git \
	curl \
    software-properties-common \
    dirmngr \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libfontconfig-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libbz2-dev \
    libblas-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && DEBIAN_FRONTEND="noninteractive" add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && DEBIAN_FRONTEND="noninteractive" apt-get update \
    && DEBIAN_FRONTEND="noninteractive" apt-get install -y --no-install-recommends r-base \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R --slave -e 'install.packages(c("tidyverse", "conflicted", "plinkQC", "RNOmni", "tclust", "gtools", "ggrepel"), repos="https://cran.rstudio.com/")'

# Download PLINK 1.9
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20240818.zip \
    && unzip -o plink_linux_x86_64_20240818.zip \
    && mv plink /usr/bin/ \
    && rm plink_linux_x86_64_20240818.zip

# Download PLINK v2.0
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241124.zip \
    && unzip -o plink2_linux_x86_64_20241124.zip \
    && mv plink2 /usr/bin/ \
    && rm plink2_linux_x86_64_20241124.zip

# Install bcftools
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git \
    && git clone https://github.com/samtools/bcftools.git \
    && cd bcftools \
    && make \
    && mv bcftools /usr/bin/ \
    && cd ..

# Install Miniconda and REGENIE
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/bin/conda \
    && rm Miniconda3-latest-Linux-x86_64.sh

RUN conda create -y -n regenie \
    && conda install -y -n regenie -c conda-forge bgenix \
    && conda install -y -n regenie -c bioconda regenie

# Add conda initialization to shell
RUN echo "source /usr/bin/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && echo "conda activate regenie" >> ~/.bashrc

SHELL ["/bin/bash", "--login", "-c"]

