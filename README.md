# EPICloneProcessing

With the Snakemake pipeline in this repository, you can process scTAM-seq data for EPI-clone.

## Installation

EPICloneProcessing requires the following software installations:

- conda, e.g., installed through [miniforge](https://github.com/conda-forge/miniforge)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- the on-premise pipeline of [Mission Bio](https://missionbio.com/). Please get in contact with your contact at Mission Bio to request access to the local pipeline and install it on your computer/server. This will create a conda environmet, which you will have to specify in the [Snakemake](Snakemake) file.

## Usage

You'll just have to fill in the right parameters in [config.yaml](configl.yaml) and then start the pipeline with:

```
snakemake --profile sge --jobs 1 --cluster-config cluster.yaml
```

This assumes that the profile sge was installed or that there is a folder called sge in the snakemake directory. For support of further compute clusters, please have a look [here](https://snakemake.readthedocs.io/en/v6.15.1/snakefiles/configuration.html).

## Contact

For questions, you can contact [Michael Scherer](mailto:michael.scherer@dkfz.de).
