[![Build Status](https://travis-ci.com/ribosomeprofiling/RiboR.svg?branch=master)](https://travis-ci.com/ribosomeprofiling/RiboR)
[![DOI](https://zenodo.org/badge/200300903.svg)](https://zenodo.org/badge/latestdoi/200300903)


# RiboR
R interface for .ribo files

## Ribo Ecosystem
The paper associated with this package and its larger ecosystem can be found [here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa028/5701654).

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

As it currently stands, this package lives on the `devel` branch, and it is not on the current release (Bioconductor 3.10). Please use the alternative installation instructions given above.

RiboR will be available via Bioconductor on the next release cycle (Bioconductor 3.11) which is scheduled for the end of April. 

Once it is available, the download instructions will be as follows.

### Install Latest Version from Bioconductor (Not Yet Available)
`if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")`

`BiocManager::install("ribor")`

## Documentation

[Here](https://ribosomeprofiling.github.io/ribor/ribor.html) is a walk-through of RiboR.

## Citing

[RiboFlow, RiboR and RiboPy: an ecosystem for analyzing ribosome profiling data at read length resolution, H. Ozadam, M. Geng, C. Cenik
Bioinformatics 36 (9), 2929-2931](https://academic.oup.com/bioinformatics/article/36/9/2929/5701654)

```bibtex
@article{ozadam2020riboflow,
  title={RiboFlow, RiboR and RiboPy: an ecosystem for analyzing ribosome profiling data at read length resolution},
  author={Ozadam, Hakan and Geng, Michael and Cenik, Can},
  journal={Bioinformatics},
  volume={36},
  number={9},
  pages={2929--2931},
  year={2020},
  publisher={Oxford University Press}
}
```
