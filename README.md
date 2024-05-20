# SH MATCHING analysis tool

[![run with singularity](https://img.shields.io/badge/run%20with-singularity-blue?style=flat&logo=singularity)](https://sylabs.io/docs/)
[![Github_Status_Badge](https://img.shields.io/badge/GitHub-2.0.0-blue.svg)](https://github.com/TU-NHM/sh_matching_pub)
[![GitHub license](https://img.shields.io/github/license/TU-NHM/sh_matching_pub)](https://github.com/TU-NHM/sh_matching_pub/blob/master/LICENSE.md)

**NB! master branch is used as development branch. Please check out [Releases](https://github.com/TU-NHM/sh_matching_pub/releases) to download a specific version of the SH matching tool.**

Developed as part of [EOSC-Nordic](https://www.eosc-nordic.eu/) project (task 5.2.1: Cross-border data processing workflows), UNITE SH matching analysis is a digital service for the global species discovery from eDNA (environmental DNA). SH matching service is based on the [UNITE](https://unite.ut.ee) datasets hosted in [PlutoF](https://plutof.ut.ee). Its output includes information about what species are present in eDNA samples, are they potentially undescribed new species, where are they found in other studies, are they alien or threatened species, etc. The output will provide DOI (Digital Object Identifier) based stable identifiers for the communicating species found in eDNA. DOIs are connected to the taxonomic backbone of [PlutoF](https://plutof.ut.ee) and [GBIF](https://www.gbif.org). In this way every DOI is accompanied by a taxon name which is still widely used for the communication of species. In the case of undescribed species, DOIs will soon be issued by the [PlutoF](https://plutof.ut.ee) system (only if SH matching service integrated with the [PlutoF](https://plutof.ut.ee) platform is used for the analysis). SH matching service covers all Eukaryota by using rDNA ITS marker sequences accompanied by sample metadata.

The script expects input files in FASTA format. Outdata files are described in [sh_matching_analysis/readme.txt](https://github.com/TU-NHM/sh_matching_pub/blob/master/sh_matching_analysis/readme.txt).

## Third-party software used by this tool

* [USEARCH](https://www.drive5.com/usearch/)
* [VSEARCH](https://github.com/torognes/vsearch)
* [ITSx](https://microbiology.se/software/itsx/)
* [KronaTools](https://github.com/marbl/Krona/wiki/KronaTools)

## Setup

### Pre-requisites

* [Singularity](https://sylabs.io/singularity/) - install Singularity (tested with version 3.5) and obtain API key for remote build

### Setup steps

1. Create Singularity Image File (SIF)
    ```console
    sudo singularity build sh_matching.sif sh_matching.def
    ```

2. OPTIONAL: Copy SIF to HPC
    ```console
    scp sh_matching.sif example_hpc_user@example.com:
    ```

3. Create input, output and working data directories
    ```console
    mkdir userdir
    mkdir indata
    mkdir outdata
    ```

4. Download FASTA dbs (https://app.plutof.ut.ee/filerepository/view/6884701) and create UDB formatted dbs
    ```console
    wget https://s3.hpc.ut.ee/plutof-public/original/d3d8b3de-83af-4fb5-b82b-359f7b730f84.zip
    mv d3d8b3de-83af-4fb5-b82b-359f7b730f84.zip sh_matching_data_udb_0_5.zip
    unzip sh_matching_data_udb_0_5.zip
    rm sh_matching_data_udb_0_5.zip
    cd data_udb/
    vsearch --makeudb_usearch sanger_refs_sh.fasta --output sanger_refs_sh.udb
    rm sanger_refs_sh.fasta
    vsearch --makeudb_usearch sanger_refs_sh_full.fasta --output sanger_refs_sh_full.udb
    rm sanger_refs_sh_full.fasta
    ```

## Running the analysis

**NB! The script expects input files in FASTA format, named as source_[run_id] and placed in indata/ directory. Outdata files are described in [sh_matching_analysis/readme.txt](https://github.com/TU-NHM/sh_matching_pub/blob/master/sh_matching_analysis/readme.txt).**

5. Run the pipeline using SIF (example data with -

* run_id=11
* region=itsfull[default]|its2
* itsx_step=yes[default]|no - flag indicating whether to include the ITSx step in the analysis (default, "yes")
* remove_userdir=yes[default]|no - flag indicating whether to delete the user directory upon pipeline completion (default, "yes")
* include_vsearch_step=yes|no[default] - flag indicating whether to include the vsearch substring dereplication step (default, "no")
* conduct_usearch_05_step=yes|no[default] - flag indicating whether to conduct the usearch complete-linkage clustering at 0.5% dissimilarity (default, "no")

    ```console
    ./sh_matching.sif /sh_matching/run_pipeline.sh 11 itsfull yes yes no no
    ```

## Citing

When using this resource, please cite as:

Abarenkov K, Kõljalg U, Nilsson RH (2022) UNITE Species Hypotheses Matching Analysis. Biodiversity Information Science and Standards 6: e93856. [https://doi.org/10.3897/biss.6.93856](https://doi.org/10.3897/biss.6.93856)

## Funding

The work is supported by [EOSC-Nordic](https://eosc-nordic.eu/) and the Estonian Research Council grant (PRG1170).
