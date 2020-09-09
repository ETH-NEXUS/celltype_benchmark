# Cell typing_benchmark
Pipeline to run individual cell typing methods and compare the accuracies of each method in a benchmark.

## Requirements:
Make sure `miniconda` is installed as `conda` environments will be used for some of the cell typing methods.
All R packages and Python modules will be installed via conda.

On Leomed:
When running each pipeline with snakemake, make sure to use the '--use-conda' argument to install and use the correct conda environment. The list of packages that will be installed in the conda environment can be found in the YML file in the `doc` folder.

## Input data format:
The data to be used as input for the benchmark needs to be an SCE file containing pre-processed counts data. Additionally, a tab-separated values (tsv) file containing a column with barcodes and a column with cell type ground-truth labels needs to be provided. Use the `add_sce_groundtruth.R` script to annotate the SCE data file with the ground-truth labels. This is the complete file to be used as input for the benchmark pipelines.

A simple sample map is used:
```
1	sample_name	N	1
```
Currently, only the sample name needs to be substituted.

## Usage:
The benchmark consists of two pipelines: a training pipeline and a testing pipeline.

For the training pipeline, the provided `training_master.snake` snake file should be used.

Once a training run has been successfully completed, you will have the input files needed to run the latter testing pipeline. For the testing pipeline, the `testing_master.snake` snake file should be used.

