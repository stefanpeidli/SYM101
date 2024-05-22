# SYMPHONY Screen Sym_101

# Running this pipeline
- Step1: install conda / mamba if necessary
- Step2: Clone this repository
- Step3: Go to `configuration/config.yaml` and change the `TEMPDIR` and `DATADIR` to the location of the data
- Step4: run `mamba env create -f environments/py_env.yaml` and activate it
- Step5: run go to the snakemake folder and run `snakemake --use-conda` to run the pipeline

# Pipeline plan for screen with single-cell readout
- Tools for single-cell readout
  - Cell Ranger
  - MAGeCK
  - Mimosca
  - BAGEL