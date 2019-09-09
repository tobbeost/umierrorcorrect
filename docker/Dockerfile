################## BASE IMAGE ######################
FROM alpine:3.7

################## METADATA ######################
LABEL base.image="alpine:3.7"
LABEL version="1"
LABEL software="umierrorcorrect"
LABEL software.version="v0.12"
LABEL about.summary="Umi error correct pipeline"
LABEL about.home="https://github.com/tobbeost/umierrorcorrect"
LABEL about.tags="Genomics"

################### Install bwa ##################
RUN apk upgrade --update && apk add --no-cache python3 python3-dev gcc gfortran freetype-dev musl-dev libpng-dev g++ lapack-dev bzip2-dev xz-dev 
#RUN apt-get install -y bwa

################## create folders ##################
USER root
RUN mkdir -p /srv/code && \
    mkdir -p /srv/example && \
    mkdir -p /srv/references

################## Install uec ##################

ADD https://github.com/tobbeost/umierrorcorrect/archive/v0.12.tar.gz /srv/code/v0.12.tar.gz
RUN pip3 install /srv/code/v0.12.tar.gz
