FROM ubuntu:16.04
MAINTAINER Haidong Yi haidyi@mail.nankai.edu.cn

# Config Software Versions here
ENV CONDA_VERSION=3-4.6.14
ENV PYTHON3_VERSION=3.6.5
ENV PYTHON2_VERSION=2.7.11

# Install some basic utilities
RUN apt-get clean && apt-get update --fix-missing && apt-get install -y \
    make \
    sudo \
    gcc \
    g++ \
    git \
    wget \
    curl \
    default-jre \
    parallel \
    cpanminus \
    bzip2 \
    libx11-6 \
    hmmer \
    emboss \
    emboss-lib \
    ncbi-blast+ \
    bioperl \
    bioperl-run \
    libdatetime-perl \
    libxml-simple-perl \
    libdigest-md5-perl \
    clustalw \
    muscle \
    prodigal \
  && rm -rf /var/lib/apt/lists/*

# Install diamond environment
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz -P /opt \
  && cd /opt && tar -xzvf diamond-linux64.tar.gz && chmod +x diamond && rm diamond-linux64.tar.gz diamond_manual.pdf \
  && mv diamond /usr/local/bin

# Install Miniconda
RUN curl -so /opt/miniconda.sh https://repo.continuum.io/miniconda/Miniconda${CONDA_VERSION}-Linux-x86_64.sh \
  && chmod +x /opt/miniconda.sh \
  && /opt/miniconda.sh -b -p /opt/miniconda \
  && rm /opt/miniconda.sh
ENV PATH=/opt/miniconda/bin:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false

# Create a python 3 and python 2 environment
RUN /opt/miniconda/bin/conda install conda-build \
  && /opt/miniconda/bin/conda create -y --name py36 python=${PYTHON3_VERSION} \
  && /opt/miniconda/bin/conda create -y --name py27 python=${PYTHON2_VERSION} \
  && /opt/miniconda/bin/conda clean -ya

# Create a soft link for python2 and pip2
ENV CONDA_DEFAULT_ENV=py36
ENV CONDA_PREFIX=/opt/miniconda/envs/${CONDA_DEFAULT_ENV}
RUN sudo ln -s /opt/miniconda/envs/py27/python ${CONDA_PREFIX}/bin/python2 \
  && sudo ln -s /opt/miniconda/envs/py27/pip ${CONDA_PREFIX}/bin/pip2
ENV PATH=${CONDA_PREFIX}/bin:${PATH}

# Install Python library acr_aca_finder needs
RUN conda install -y biopython && conda clean -ya

# set app directory
RUN mkdir -p /app

# Install the acr_aca_finder and CRISPRCas-Finder
RUN cd /app && git clone https://haidyi:yhd19930426@github.com/haidyi/acr_aca_finder.git
RUN cd /app/acr_aca_finder/dependencies/CRISPRCasFinder/ && chmod +x installer_UBUNTU.sh && ./installer_UBUNTU.sh

# Config some environmental varialbes for CRISPRCas-Finder
RUN sed -i '1c #!/usr/bin/env python2' /app/acr_aca_finder/dependencies/CRISPRCasFinder/macsyfinder-1.0.5/bin/macsyfinder \
  && sed -i '$c export PATH=/app/acr_aca_finder/dependencies/CRISPRCasFinder/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH' ~/.profile
ENV MACSY_HOME=/app/acr_aca_finder/dependencies/CRISPRCasFinder/macsyfinder-1.0.5/
ENV PATH=/app/acr_aca_finder/dependencies/CRISPRCasFinder/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# Set the working directory
WORKDIR /app/acr_aca_finder

# CMD
CMD [ "python3 acr_aca_cri_runner.py -n sample_organisms/GCF_000210795.2/GCF_000210795.2_genomic.fna -o test_output -c 0 -z B"]