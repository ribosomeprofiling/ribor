[![Build Status](https://travis-ci.com/ribosomeprofiling/RiboR.svg?branch=master)](https://travis-ci.com/ribosomeprofiling/RiboR)
[![DOI](https://zenodo.org/badge/200300903.svg)](https://zenodo.org/badge/latestdoi/200300903)


# RiboR
R interface for .ribo files

## Installation

RiboR requires R version **3.6** or higher.



### Note For Linux Users

Install dependencies for devtools.
For Ubuntu based distributions, you can use the following command.

`sudo apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev build-essential m4 autoconf -y`

### Install Latest Version From GitHub

1) Install devtools

`install.packages("devtools")`

2) Load the devtools package

`library("devtools")`

3) Install RiboR from github

`install_github("ribosomeprofiling/ribor")`

### Note for Bioconductor Installation

For now, RiboR is **NOT** a bioconductor package. Please use the alternative installation instructions given above.

RiboR has been submitted to Bioconductor for evaluation. We hope that RiboR will be available via Bioconductor on the next release cycle.
Once it is available, the download instructions will be as follows.

### Install Latest Version from Bioconductor (Not Yet Available)
`if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")`

`BiocManager::install("ribor")`

## Documentation

[Here](https://ribosomeprofiling.github.io/ribor/ribor.html) is a walk-through of RiboR.
