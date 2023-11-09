# EPICloneProcessing
With the Snakemake pipeline in this repository, you can process scTAM-seq data for EPI-clone.

## Usage
After installing the required software packages (the pipeline should automatically run on the CRG cluster, given that no configuration changed), you'll just have to fill in the right parameters in [config.yaml](configl.yaml) and then start the pipeline with:

```
snakemake --profile sge --jobs 1 --cluster-config cluster.yaml
```

This assumes that the profile sge was installed or that there is a folder called sge in the snakemake directory.
