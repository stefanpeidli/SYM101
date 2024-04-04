# SYMPHONY

TODO: check for collateral activity in Wessels et al data!

# Running this pipeline
- Step1: install conda / mamba if necessary
- Step2: Clone this repository
- Step3: Go to `configuration/config.yaml` and change the `TEMPDIR` and `DATADIR` to the location of the data
- Step4: run `mamba env create -f environments/py_env.yaml` and activate it
- Step5: run go to the snakemake folder and run `snakemake --use-conda` to run the pipeline

# Plan for eval and QC (From OP)
- Check mitochondrial reads targeted vs NT. Should be the same roughly.
- Average strength of transcriptomic knockdown reported to be mean 65% with s.d. 8.7%

# Pipeline plan for screen with single-cell readout
- Tools for single-cell readout
  - Cell Ranger
  - MAGeCK
  - Mimosca
  - BAGEL