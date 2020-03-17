################## BASE IMAGE ######################
FROM ubuntu:18.04

################## METADATA ######################
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="umierrorcorrect"
LABEL software.version="v0.15"
LABEL about.summary="Umi error correct pipeline"
LABEL about.home="https://github.com/tobbeost/umierrorcorrect"
LABEL about.tags="Genomics"

# Add umierrorcorrect download from github to container in /etc
ADD https://github.com/tobbeost/umierrorcorrect/archive/v0.15.tar.gz /etc/

# Update Ubuntu Software repository and install dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    make \
    curl \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    bzip2 \
    xz-utils \
    gcc \
    g++ \
    python3-dev \
    python3-pip \
    python3-setuptools \
    pigz \
    #samtools \
    bwa && \
    pip3 install wheel && \
    pip3 install pysam && \
    pip3 install scipy && \
    pip3 install matplotlib && \
    pip3 install /etc/v0.15.tar.gz


################## create folders ##################
#USER root
#RUN mkdir -p /srv/code && \
#    mkdir -p /srv/example && \
#    mkdir -p /srv/references
#
################## Install uec ##################

#ADD https://github.com/tobbeost/umierrorcorrect/archive/v0.12.tar.gz /srv/code/v0.12.tar.gz
#RUN pip3 install /srv/code/v0.13.tar.gz
