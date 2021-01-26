# sh_matching
EOSC-Nordic SH matching service preparations

## Setup

### Pre-requisites

* [Singularity](https://sylabs.io/singularity/) - install Singularity and obtain API key for remote build
* [USEARCH](https://) - download USEARCH binary (usearch11.0.667_i86linux32) and place it in the main folder

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
    mkdir userdir
    mkdir indata
    mkdir outdata
    ```

## SH matching: running the analysis (OPTIONAL: This could be done through sbatch slurm scripts)

**NB! The script expects input files in FASTA format, named as source_[run_id] and placed in indata/ directory.**

4. Run the pipeline using SIF (example data with run_id=11)
    ```console
    ./sh_matching.sif /sh_matching/run_pipeline.sh 11
    ```
