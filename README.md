# ADNI sex differences analysis

This is the repository for the sex differences analysis on the ADNI dataset.
The purpose is to detect sex differences in metabolite-phenotype associations.

## Clone repository

In your terminal type

```bash
git clone https://github.com/tomszar/ADNI_project.git
```

## Environment setup

The repository provides an `environment.yml` file to use with conda

First, install [anaconda](https://www.anaconda.com/products/individual) on your local computer or user server account following their directions, using the corresponding installer.
Next, on the terminal, install the conda environment by running:

```bash
conda env create -f environment.yml
```

## Replicate the analysis

Depending on whether you're replicating the analysis on your local machine or sending the job to a cluster, you can run either `run_local.sh` or `run_cluster.sh`

### On cluster

On your server terminal, run the command

```bash
qsub run_cluster.sh
```

The parameters used in the file are the ones used in the Penn State Roar server.
Depending on the system you will need to modify, remove, or add parameters to the file.

### On local machine

Make sure the environment is activated before running by typing:

```bash
conda activate adni_project
```

Then, type:

```bash
bash run_local.sh
```
