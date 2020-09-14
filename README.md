# Cell typing_benchmark
Pipeline to run individual cell typing methods and compare the accuracies of each method in a benchmark.

## Requirements:
Make sure `miniconda` is installed as `conda` environments will be used for some of the cell typing methods.
All R packages and Python modules will be installed via conda automatically using the pre-provided YAML file.

When running each pipeline with snakemake, make sure to use the `--use-conda` argument to install and use the correct conda environment. The list of packages that will be installed in the conda environment can be found in the YML file in the `doc` folder.



Alternatively, if you do not wish to use conda, the following R and Python packages need to be installed into your working environment:

#### R packages (R version 3.6.1)

```
optparse
SingleCellExperiment
tidyverse
e1071
Hmisc
cowplot
caret
monocle3
garnett
limma
org.Hs.eg.db
Hmisc
igraph
pheatmap
randomForest
RColorBrewer
reshape2
scran
Seurat
uwot
```

#### Python modules (Python 3.7.7)

```
N/A (Base Python only)
```

**TensorFlow** also needs to be installed for the fcNN method.



##### On Leomed:
See **Usage** below for an example of how this would be run.

In the future, a pre-installed conda environment will be provided on Leomed so that a new conda environment will not be installed for each new run of the pipelines.

In the meanwhile, a conda environment can be used as follows:
`conda activate /cluster/work/nexus/marcus/miniconda3/envs/garnett/`
However, using this environment means you must not use the `--use-conda` mode of Snakemake.

## Input data:
Required input data:

- RDS file containing a SCE (SingleCellExperiment) object with pre-processed counts data & per-cell major and minor cell type (ground-truth) labels.

  

**If this is your first time running the cell type benchmark, it is unlikely that your counts data object will be annotated with correctly formatted cell type labels.** If so, please follow the steps below to first prepare your data:

### Generating a SCE counts object with ground-truth cell type labels

Required input data:

- RDS file containing a SCE (SingleCellExperiment) object with pre-processed counts data
- A *tsv* (tab-separated values) file containing a column with barcodes and a column with cell type ground-truth labels

The ground-truth *tsv* file should be formatted as follows:

- 4 columns:
  1. Barcodes (in data)
  2. Cell types
  3. Standardised major cell types
  4. Standardised minor cell types

```
Barcode	Original_type	Coerced_major_type	Coerced_minor_type
0000e539a7dfe057f8013c9f5081b369	naive B cell	B.cells	B.cells.naive
0002078f929f37ffa735bd59379e35f7	T-helper cell	T.cells	T.cells.CD4
00032d6cd1f7c21a3ac0ffb1b0233a61	precursor B cell	B.cells	B.cells.precursor
00050b3d3f47102bd42ee12caf0f7320	cytotoxic T cell type 2	T.cells	T.cells.CD8
00065affcfb2b5b56902a4ce2cfa21aa	pro-B cell	B.cells	B.cells.precursor
000d2b88671106c8f71ab4f1586378fa	natural killer cell	NK.cells	none
0025a6ca5ecb85790da909c564508a23	T-helper cell	T.cells	T.cells.CD4
002d84f0415656febbbfaf17239cde7c	naive B cell	B.cells	B.cells.naive
002e791afa847df6290736b88435bcaf	cytotoxic T cell type 2	T.cells	T.cells.CD8
```

*Note: While the labels in the header are unimportant, a header is necessary.*

The standardised cell type column should consist of common, mutual labels that can be used to compare results between different samples. For single-sample analyses, the same number of columns are required, but they can be arbitrary as they need not match another sample's labels.



The process of generating the list above can be automated using the `generate_ground-truth_list.py` Python script.

**Usage**: `python generate_ground-truth_list.py --original_labels metadata.txt --label_dictionary x --subtype_dictionary major_subtype_mapping.txt --output formatted_ground_truth.tsv `

Example of `--original_labels` input file:

- A 2-column tsv file consisting of a barcode column and a cell type label column.

```
0000e539a7dfe057f8013c9f5081b369	naive B cell
00032d6cd1f7c21a3ac0ffb1b0233a61	precursor B cell
000a156d3ff1f1f13f30601a6dc66e72	CD14+ monocyte type 1
```

Example of `--label_dictionary` input file:

```
Dictionary_original     Dictionary_matched
cd14_monocytes_fastqs   Monocytes
CD14+ monocyte type 1   Monocytes
CD14+ monocyte type 2   Monocytes
CD16+ monocyte  Monocytes
cd4_t_helper_fastqs     T.cells.CD4
```

- A 2-column tsv file with the original cell type labels and standardized cell type labels. Like stated above, if this is not necessary, the columns can be identical.

Example of `--subtype_dictionary` input file:

- A file mapping the relationship between major cell types and minor subtypes.

An example can be found in `marker_files/PBMC/celltype_config_PBMC.txt`. 

```
Major_type      Subtypes
Monocytes       none
T.cells T.cells.CD4,T.cells.CD8,T.cells.regulatory
Dendritic.cells none
B.cells B.cells.memory,B.cells.naive,B.cells.precursor
Plasma.cells    none
Plasmacytoid.dendritic.cells    none
NK.cells        none
Stem.cells      none
```



Once the final ground-truth file has been generated, use the `add_sce_groundtruth.R` script to annotate the SCE data file with the ground-truth labels.

**Usage:** `Rscript add_sce_groundtruth.R --sce_in processed_sce.RDS --labels_in ground_truth.tsv --matched_labels --sce_out final_annotated_sce.RDS ` 

This output of this is the complete file to be used as input for the benchmark pipelines.

---

A simple sample map is used:
```
1       sample_name     N       1
```
Currently, only the sample name needs to be substituted.

---

Example configuration files for each Snakemake pipeline can be found in `doc`.

- `config_training_pipeline.yml`
- `config_testing_pipeline.yml`

---

## Usage:
The benchmark consists of two pipelines: a training pipeline and a testing pipeline.

For the training pipeline, the provided `training_master.snake` snake file should be used.

Once a training run has been successfully completed, you will have the input files needed to run the latter testing pipeline. For the testing pipeline, the `testing_master.snake` snake file should be used.

**Leomed:**

Modules to be loaded first:

- `module load python_cpu/3.7.4 `

Snakemake installation that will be used:

- `/cluster/apps/python/3.7.4_cpu/x86_64/bin/snakemake`

Example dry run:

- ```sh
  snakemake --notemp --latency-wait 36000 -s git/snake/training_master.snake --configfile git/doc/config_training_pipeline.yml --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k -n
  ```

Example training run:

- ```sh
  bsub -J ct_benchmark_train -W 23:59 -R "rusage[mem=2000]" -eo ct_benchmark_train.err -oo ct_benchmark_train.out "snakemake --notemp --latency-wait 36000 -s git/snake/training_master.snake --configfile git/doc/config_training_pipeline.yml --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k"
  ```

Example testing run:

- ```sh
  bsub -J ct_benchmark_test -W 23:59 -R "rusage[mem=2000]" -eo ct_benchmark_test.err -oo ct_benchmark_test.out "snakemake --notemp --latency-wait 36000 -s git/snake/testing_master.snake --configfile git/doc/config_testing_pipeline.yml --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k"
  ```

  
