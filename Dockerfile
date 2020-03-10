FROM ubuntu:18.04

# Config Software Versions here
ENV CONDA_VERSION=3-4.6.14
ENV PYTHON3_VERSION=3.6.5

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
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.26/diamond-linux64.tar.gz -P /opt \
  && cd /opt && tar -xzvf diamond-linux64.tar.gz && chmod +x diamond && rm diamond-linux64.tar.gz diamond_manual.pdf \
  && mv diamond /usr/local/bin

# Install Miniconda
RUN curl -so /opt/miniconda.sh https://repo.continuum.io/miniconda/Miniconda${CONDA_VERSION}-Linux-x86_64.sh \
  && chmod +x /opt/miniconda.sh \
  && /opt/miniconda.sh -b -p /opt/miniconda \
  && rm /opt/miniconda.sh
ENV PATH=/opt/miniconda/bin:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false

# Create a python 3 environment
RUN /opt/miniconda/bin/conda install conda-build \
  && /opt/miniconda/bin/conda create -y --name py36 python=${PYTHON3_VERSION} \
  && /opt/miniconda/bin/conda clean -ya

# Config the environment variable for python3
ENV CONDA_DEFAULT_ENV=py36
ENV CONDA_PREFIX=/opt/miniconda/envs/${CONDA_DEFAULT_ENV}
ENV PATH=${CONDA_PREFIX}/bin:${PATH}

# Install Python library - biopython acr_aca_finder needs
RUN conda install -y biopython && conda clean -ya

# set app directory
RUN mkdir -p /app

# Install the acr_aca_finder and CRISPRCas-Finder
RUN cd /app && git clone https://github.com/haidyi/acrfinder.git
RUN cd /app/acrfinder/dependencies/CRISPRCasFinder/ && chmod +x installer_UBUNTU.sh && ./installer_UBUNTU.sh

# make prophage database
RUN cd /app/acrfinder/dependencies/prophage && makeblastdb -in prophage_virus.db -dbtype prot -out prophage

# make cdd database
RUN mkdir -p /app/acrfinder/dependencies/cdd
RUN cd /app/acrfinder/dependencies/cdd && wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz && tar -xzf cdd.tar.gz && rm cdd.tar.gz
RUN cd /app/acrfinder/dependencies/cdd && makeprofiledb -title CDD.v.3.12 -in Cdd.pn -out Cdd -threshold 9.82 -scale 100.0 -dbtype rps -index true

# make cdd-mge database
RUN cd /app/acrfinder/dependencies/ && tar -xzf cdd-mge.tar.gz && rm cdd-mge.tar.gz

# Config some environmental varialbes for CRISPRCas-Finder
RUN sed -i '1c #!/usr/bin/env python2' /app/acrfinder/dependencies/CRISPRCasFinder/macsyfinder-1.0.5/bin/macsyfinder \
  && sed -i '$c export PATH=/app/acrfinder/dependencies/CRISPRCasFinder/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH' ~/.profile
ENV MACSY_HOME=/app/acrfinder/dependencies/CRISPRCasFinder/macsyfinder-1.0.5/
ENV PATH=/app/acrfinder/dependencies/CRISPRCasFinder/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# Set the working directory
WORKDIR /app/acrfinder

# CMD
CMD ["python3 acr_aca_cri_runner.py -n sample_organisms/GCF_000210795.2/GCF_000210795.2_genomic.fna -o test_output -c 0 -z B -c 2 -p true -g true"]