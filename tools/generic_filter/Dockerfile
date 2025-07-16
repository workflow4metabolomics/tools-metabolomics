################################################################################
### 
### [CONTAINER CORE FUNCTIONS]: 
###     install "Tool - Generic Filter" Galaxy tool (and required third part softwares, libraries, ...).
### [NOTE]
###     please refer to README.md and about_docker.md files for further informations
### 
################################################################################

################################################################################
### fix parent containter
FROM ubuntu:16.04

################################################################################
### set author
MAINTAINER Nils Paulhe <nils.paulhe@inra.fr>

################################################################################
### sets the environment variables
ENV TOOL_VERSION = "release_2016.03.03"
ENV CONTAINER_VERSION = 0.1

LABEL version = "${CONTAINER_VERSION}"
LABEL tool_version = "${TOOL_VERSION}"

################################################################################
### install third part tools 

# add debian repo for latest version of R
RUN echo "deb http://cran.univ-paris1.fr/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 

# Update and upgrade system
RUN apt-get update && \
    apt-get -y upgrade

# install R
RUN apt-get install -y \ 
    r-base \
    libcurl4-openssl-dev \
    libxml2-dev
# NOTE: add `apt-get install -y git` if required

# init R env. (Docker)
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

# install R libs
RUN Rscript -e "install.packages('batch', dep=TRUE)"

################################################################################
### install core scripts

# init. WORKDIR
RUN [ "mkdir", "/scripts" ]

#
# [NOTE] to add scripts, we have two options: get them from GitHub OR copy them from this directory
# 

# get scripts using Git (option 1)
# RUN cd /scripts && \
#    git clone -b release/${TOOL_VERSION} --recursive https://github.com/workflow4metabolomics/tool-generic_filter.git

# copy scripts files from this directory (option 2)
COPY "." "/scripts/"

## set WORKDIR
# WORKDIR "/scripts"

# set authorizations
RUN ["chmod", "a+x", "/scripts/filter_wrap.R"]

# make tool accessible through PATH
ENV PATH = $PATH:/scripts

################################################################################
### clean
RUN apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/ /tmp/* /var/tmp/*
# NOTE: run `apt-get remove -y git && \` if required

################################################################################
### Define Entry point script
## ENTRYPOINT ["/scripts/filter_wrap.R"]

### [END]