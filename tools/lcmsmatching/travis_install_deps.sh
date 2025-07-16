#!/bin/bash

# Install devtools
R -e "install.packages('devtools', dependencies = TRUE, repos = 'https://cloud.r-project.org/')"

# Get all dependencies except biodb
dependencies=$(grep 'requirement.*type="package.*>' lcmsmatching.xml | grep -v biodb | sed 's!^.*>r-\([^<>]*\)</.*$!\1!')

# Install dependencies
for dep in $dependencies ; do

	# Get version
	version=$(grep "requirement.*type=\"package.*>r-$dep.*R_VERSION=" lcmsmatching.xml | sed 's/^.*R_VERSION="\([^"]*\)".*$/\1/')
	if [[ -z $version ]] ; then
		version=$(grep "requirement.*type=\"package.*>r-$dep" lcmsmatching.xml | sed 's/^.*version="\([^"]*\)".*$/\1/')
	fi

	# Install package
	R -e "devtools::install_version('$dep', version = '$version', repos = 'https://cloud.r-project.org/')"
done

# Install right version of biodb
version=$(grep r-biodb lcmsmatching.xml | sed 's/^.*version="\([^"]*\)".*$/\1/')
R -e "devtools::install_github('pkrog/biodb', ref = 'v$version')"
