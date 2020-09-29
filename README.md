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

  

**If this is your first time running the cell type benchmark, it is unlikely that your counts data object will be annotated with correctly formatted cell type labels.**

See the "Data" section near the end of the Readme for examples of entire workflows in preparing data (such as the example data used) for use in the benchmark.

If your data is already processed, please follow the steps below to first prepare your data:

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

# Example data - processing workflows

There are 3 data sets provided for an example run of the cell typing benchmark.
 
- Zheng (sorted & merged)
- (2020 Broad) Adult PBMC
- (2020 Broad) Newborn PBMC

Their pre-processing, processing, and annotation steps are detailed below.

Unless otherwise specified, the paths link to files on Leomed.

##########################################
            Table of contents            
##########################################
[1] Zheng
    (a) Pre-processing (cellranger count & aggr)
    (b) Processing
    (c) Subset
        i. Processing
    (d) Ground-truth
    (e) Cell type proportions table
    (f) Summary of file paths
[2] Adult PBMC
    (a) Processing
    (b) Subset
        i. Processing
    (c) Ground-truth
    (d) Cell type proportions table
    (e) Summary of file paths
[3] Newborn PBMC
    (a) Processing
    (b) Subset
        i. Processing
    (c) Ground-truth
    (d) Cell type proportions table
    (e) Summary of file paths

##########################################
Pre-processing and processing necessary for each respective datase
##########################################

## [1] Zheng data set
Consists of 10 FACS sorted immune cell populations which are then merged into one dataset.
Data sets:
1. CD19+ B cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/b_cells_fastqs/fastqs
2. CD14+ Monocytes
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/cd14_monocytes_fastqs/fastqs
3. CD34+ Cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/cd34_fastqs/fastqs
4. CD4+ Helper T cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/cd4_t_helper_fastqs/fastqs
5. CD56+ Natural Killer Cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/cd56_nk_fastqs/fastqs
6. CD8+ Cytotoxic T cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/cytotoxic_t_fastqs/fastqs
7. CD4+/CD45RO+ Memory T cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/memory_t_fastqs/fastqs
8. CD8+/CD45RA+ Naive Cytotoxic T cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/naive_cytotoxic_fastqs/fastqs
9. CD4+/CD45RA+/CD25- Naive T cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/naive_t_fastqs/fastqs
10. CD4+/CD25+ Regulatory T cells
Path to FASTQs: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/regulatory_t_fastqs/fastqs

FASTQ files for each dataset were retrieved from the '10X Genomics Single Cell Gene Expression Datasets' page.
(https://support.10xgenomics.com/single-cell-gene-expression/datasets - Single Cell 3' Paper: Zheng et al. 2017)

----------------------------------------------------------------------------------------------------------------
### (a) Pre-processing - cellranger count

Each dataset was processed separately with the cellranger count pipeline.
Example of commands used (on Euler) - one per dataset above:
```
module load /cluster/project/nexus/utilities/sharedPrograms/cellranger/cellranger-3.1.0/cellranger.modulefile;
cellranger count --id=b_cells --transcriptome=/cluster/project/nexus/utilities/databases/singlecell/10xGenomics/gene_expression_3_0_2/refdata-cellranger-GRCh38-3.0.0 --localcores=24 --fastqs=/cluster/project/nexus/benchmarking/celltyping/data/zheng_sorted/fastqs/b_cells_fastqs/fastqs/flowcell1/ --nosecondary
```
Output location: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/Zheng_datasets_cellranger

----------------------------------------------------------------------------------------------------------------
#### Merging 10 datasets - cellranger aggr

Using the output of the cellranger count pipeline from each of the 10 datasets, cellranger aggr was used to merge the data into a single dataset.

To avoid issues of duplicated barcodes, the barcode of each respective dataset is appended with a unique suffix using the numbering in the list above. (E.g. All B cells barcodes will have a '-1' added to the end, '-2' for monocytes, and so on.)
This can be done through a simple command as follows:
`sed -i s/$/-1/ barcodes.tsv`
(This will append '-1' to the end of each line/barcode.)

A csv file has to be provided as input mapping each population to its 'molecule_info.h5' file from its respective cellranger count output, meaning one line per dataset.
Example showing header and line for one dataset:
```
library_id,molecule_h5,FACS_type
b_cells,/cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/Zheng_datasets_cellranger/b_cells/outs/molecule_info.h5,b_cells
```

Command used (on Euler): `cellranger aggr --id=zheng_merged --csv=aggregate_libraries.csv --normalize=none`
(NOTE: Command was run from within output directory below.)
Output location: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/Zheng_cellranger_merge/cellranger_aggr_output/outs

----------------------------------------------------------------------------------------------------------------
### (b) Processing - running scRNAseq pipeline on dataset

Snakemake command used (on Euler):
```
export HDF5_USE_FILE_LOCKING=FALSE; /cluster/project/nexus/utilities/sharedPrograms/snakemake/snakemake_4.8.0/bin/snakemake --notemp --latency-wait 240 -s /cluster/project/nexus/marcus/git/singlecell_analysis/snake_master/snake_gex_master.snake --configfile /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/snake_analysis_files/config.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R rusage[mem={params.mem},scratch={params.scratch}] -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k;
```

Output location: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/analysis

Prior to the run, barcodes.tsv, matrix.mtx, and features.tsv files taken from the cellranger aggr run output in /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/Zheng_cellranger_merge/cellranger_aggr_output/outs/filtered_feature_bc_matrix were copied to /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/analysis/cellranger_run_gex/.

##########################################

### (c) Generating subset

The first step of the scRNAseq processing pipeline is generating an hdf5 file from the cellranger count output containing the raw count data.
This is used to generate an hdf5 file containing a subset of the original data. Around 10% of cells in the original dataset are chosen randomly, respecting and maintaing the proportion of cell types in the dataset.

Subsetting was performed with an R script as follows (on Euler):
`Rscript /cluster/project/nexus/marcus/git/celltype_benchmark/scripts/subset_data_proportional.R --input_data /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/analysis/rawCounts/Zheng_merged.h5 --true_labels /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/subset/zheng_ground_truth.tsv --output_subset /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_hdf5/Zheng_subset.h5`

Output location: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_hdf5/

----------------------------------------------------------------------------------------------------------------
#### i. Processing - running scRNAseq pipeline on subset

Snakemake command used (on Euler):
```
export HDF5_USE_FILE_LOCKING=FALSE; /cluster/project/nexus/utilities/sharedPrograms/snakemake/snakemake_4.8.0/bin/snakemake --notemp --latency-wait 240 -s /cluster/project/nexus/marcus/git/singlecell_analysis/snake_master/snake_gex_master.snake --configfile /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/snake_analysis_files/config.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R rusage[mem={params.mem},scratch={params.scratch}] -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k;
```

Output location: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_scrun/analysis/

Prior to the run, a symbolic link to the subset hdf5 file is created in the pipeline run directory:
`ln -s /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_hdf5/Zheng_subset.h5 /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_scrun/analysis/rawCounts/Zheng_subset.h5`

##########################################

### (d) Annotating data with cell type ground-truth labels

For the benchmark, the data must first be labeled with the true cell types per barcode.
To do this, a file (hereafter referred to as a ground-truth dictionary) mapping the labels to each barcode must be generated.
Note that this only needs to be done once per (full) data set.
Please see the documentation on the celltyping_benchmark Git repo for examples and required formatting of the input files needed.

Command used to generate ground-truth dictionary for the Zheng data:
`python /cluster/work/nexus/benchmarking/celltyping/git/scripts/generate_ground-truth_list.py --original_labels /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/Zheng_cellranger_merge/cellranger_aggr_output/outs/filtered_feature_bc_matrix/barcodes.tsv --label_dictionary /cluster/work/nexus/benchmarking/celltyping/git/doc/matched_types_dictionary.tsv --output /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/zheng_ground_truth.tsv --subtype_dictionary /cluster/work/nexus/benchmarking/celltyping/data/major_subtype_mapping.txt`

----------------------------------------------------------------------------------------------------------------

Once this is done, the processed RDS file from the scRNAseq processing pipeline's 'prep_celltyping' step is ready to be annotated.

Command used to annotate processed RDS with ground-truth labels (Zheng data):
```
Rscript /cluster/work/nexus/benchmarking/celltyping/git/scripts/add_sce_groundtruth.R --sce_in /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/analysis/prep_celltyping/Zheng_merged.genes_cells_filtered.corrected.RDS --sce_out /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/extended_sce/Zheng_merged_annotated.RDS --labels_in /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/zheng_ground_truth.tsv
```

Output path: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/extended_sce/

----------------------------------------------------------------------------------------------------------------

Repeat the annotation process, but this time using the subset data.

Command used to annotate processed RDS with ground-truth labels (Zheng data):
```
Rscript /cluster/work/nexus/benchmarking/celltyping/git/scripts/add_sce_groundtruth.R --sce_in /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_scrun/analysis/prep_celltyping/Zheng_subset.genes_cells_filtered.corrected.RDS --sce_out /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/Zheng_subset_annotated.RDS --labels_in /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/zheng_ground_truth.tsv
```

Output path: /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/

##########################################

### (e) Table of cell type proportions in dataset (full & subset)

A script was used to generated tables summarizing each dataset's cell type counts and proportions.
The prefixes of the column labels map as follows:
- "raw_full": barcodes file (cellranger output)
- "raw_subset": HDF5 file of subset data (file generated by subsetting script)
- "processed_full": RDS file of SCE object (output of scRNAseq pipeline: 'prep_celltyping')
- "processed_subset": RDS file of subsetted SCE object (output of scRNAseq pipeline: 'prep_celltyping')
- "metadata": ground-truth dictionary

Command used to generated Zheng summary table (on Leomed):
```
Rscript /cluster/work/nexus/benchmarking/celltyping/git/scripts/write_celltype_proportion_tables.R --data_full /cluster/dataset/nexus/benchmarking/celltyping/data/zheng_sorted/extended_sce/Zheng_merged_annotated.RDS --data_raw /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/Zheng_cellranger_merge/cellranger_aggr_output/outs/filtered_feature_bc_matrix/barcodes.tsv --data_subset /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/Zheng_subset_annotated.RDS --data_subset_raw /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/zheng_subset_barcodes.tsv --data_dict /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/zheng_ground_truth.tsv --output /cluster/dataset/nexus/benchmarking/celltyping/data/summary_tables/Zheng/ --sample_name Zheng
```

Output directory: /cluster/work/nexus/benchmarking/celltyping/data/summary_tables

##########################################

### (f) Summary of paths to relevant Zheng files

- Barcodes: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/Zheng_cellranger_merge/cellranger_aggr_output/outs/filtered_feature_bc_matrix/barcodes.tsv
- Ground-truth: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/zheng_ground_truth.tsv

- Processing run directory: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/
- Processing run (subset) directory: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_scrun/

- Data hdf5 (raw counts): /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/analysis/rawCounts/Zheng_merged.h5
- Subset hdf5 (raw counts): /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_hdf5/Zheng_subset.h5

- Processed RDS: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_merged_scrun/analysis/prep_celltyping/Zheng_merged.genes_cells_filtered.corrected.RDS
- Processed RDS (subset): /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_scrun/analysis/prep_celltyping/Zheng_subset.genes_cells_filtered.corrected.RDS

- Annotated RDS: /cluster/work/nexus/benchmarking/celltyping/data/zheng_sorted/extended_sce/Zheng_merged_annotated.RDS
- Annotated RDS (subset): /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/Zheng_subset_annotated.RDS

- Proportions table: /cluster/work/nexus/benchmarking/celltyping/data/summary_tables/Zheng/Zheng_true_labels_table.tx

----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------

## [2] 2020-Mar-Census-Adult-Immune-10x datase
Adult PBMC dataset (preprocessed) retrieved from the Broad Institute.
See /cluster/work/nexus/benchmarking/celltyping/data/broad/readme.txt for more details.

Since the data is already preprocessed (cellranger output), we can proceed directly to the processing step.

----------------------------------------------------------------------------------------------------------------
### (a) Processing - running scRNAseq pipeline on dataset

WARNING: full dataset requires a special implementation of the processing pipeline due to the size of the dataset.
As such, only the subset will be used for annotation.

Commands are still experimental and all files can be found on Euler.

Snakemake command used (on Euler):
```
export HDF5_USE_FILE_LOCKING=FALSE; /cluster/project/nexus/utilities/sharedPrograms/snakemake/snakemake_4.8.0/bin/snakemake --notemp --latency-wait 240 -s /cluster/project/nexus/benchmarking/celltyping/test_pilot/2020-Mar-Census-Adult-Immune-10x_run/snake_analysis_files/snake_gex_master.snake --configfile /cluster/project/nexus/benchmarking/celltyping/test_pilot/2020-Mar-Census-Adult-Immune-10x_run/snake_analysis_files/config.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R rusage[mem={params.mem},scratch={params.scratch}] -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k;
```

Output location: /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_run

##########################################

### (b) Generating subset

The first step of the scRNAseq processing pipeline is generating an hdf5 file from the cellranger count output containing the raw count data.
This is used to generate an hdf5 file containing a subset of the original data. Around 10% of cells in the original dataset are chosen randomly, respecting and maintaing the proportion of cell types in the dataset.

Subsetting was performed with an R script as follows (on Euler):
`Rscript /cluster/project/nexus/marcus/git/celltype_benchmark/scripts/subset_data_proportional.R --input_data /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_run/analysis/rawCounts/2020-Mar-Census-Adult-Immune-10x.h5 --true_labels /cluster/project/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x.ground_truth.tsv --output_subset /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/2020-Mar-Census-Adult-Immune-10x.scp_subset.h5`

Output location: /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/

----------------------------------------------------------------------------------------------------------------
## i. Processing - running scRNAseq pipeline on subset

Snakemake command used (on Euler):
```
export HDF5_USE_FILE_LOCKING=FALSE; /cluster/project/nexus/utilities/sharedPrograms/snakemake/snakemake_4.8.0/bin/snakemake --notemp --latency-wait 240 -s /cluster/project/nexus/benchmarking/celltyping/test_pilot/2020-Mar-Census-Adult-Immune-10x_run/snake_analysis_files/snake_gex_master.snake --configfile /cluster/project/nexus/benchmarking/celltyping/test_pilot/adult_subset/snake_analysis_files/config.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R rusage[mem={params.mem},scratch={params.scratch}] -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k;
```

Output location: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_subset_scrun

Prior to the run, a symbolic link to the subset hdf5 file is created in the pipeline run directory:
`ln -s /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/2020-Mar-Census-Adult-Immune-10x.scp_subset.h5 /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_subset_scrun/analysis/rawCounts/2020-Mar-Census-Adult-Immune-10x_subset.scp.h5`

##########################################

## (c) Annotating data with cell type ground-truth labels

For the benchmark, the data must first be labeled with the true cell types per barcode.
To do this, a file (hereafter referred to as a ground-truth dictionary) mapping the labels to each barcode must be generated.
Note that this only needs to be done once per (full) data set.
Please see the documentation on the celltyping_benchmark Git repo for examples and required formatting of the input files needed.

Command used to generate ground-truth dictionary for the Adult PBMC data:
`python /cluster/work/nexus/benchmarking/celltyping/git/scripts/generate_ground-truth_list.py --original_labels /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x_annotated_v1.scp.metadata.txt --label_dictionary /cluster/work/nexus/benchmarking/celltyping/git/doc/matched_types_dictionary.tsv --output /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x.ground_truth.tsv --subtype_dictionary /cluster/work/nexus/benchmarking/celltyping/data/major_subtype_mapping.txt`

----------------------------------------------------------------------------------------------------------------

Once this is done, the processed RDS file from the scRNAseq processing pipeline's 'prep_celltyping' step is ready to be annotated.

Command used to annotate processed RDS with ground-truth labels (Adult PBMC subset):
```
Rscript /cluster/work/nexus/benchmarking/celltyping/git/scripts/add_sce_groundtruth.R --sce_in /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_subset_scrun/analysis/prep_celltyping/2020-Mar-Census-Adult-Immune-10x_subset.scp.genes_cells_filtered.corrected.RDS --sce_out /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/2020-Mar-Census-Adult-Immune-10x_subset_annotated.RDS --labels_in /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x_annotated_v1.scp.metadata.tx
```

Output path: /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/

##########################################

## (d) Table of cell type proportions in dataset (full & subset)

A script was used to generated tables summarizing each dataset's cell type counts and proportions.
The prefixes of the column labels map as follows:
- "raw_full": barcodes file (cellranger output)
- "raw_subset": HDF5 file of subset data (file generated by subsetting script)
- "processed_subset": RDS file of subsetted SCE object (output of scRNAseq pipeline: 'prep_celltyping')
- "metadata": ground-truth dictionary
("processed_full" omitted due to size of file/object.)

Command used to generated Adult PBMC summary table (on Leomed):
```
Rscript /cluster/work/nexus/benchmarking/celltyping/data/summary_tables/write_celltype_proportion_tables.R --data_raw /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_scrun/analysis/rawCounts/2020-Mar-Census-Adult-Immune-10x.h5 --data_subset /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/2020-Mar-Census-Adult-Immune-10x_subset_annotated.RDS --data_subset_raw /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/2020-Mar-Census-Adult-Immune-10x.scp_subset.h5 --data_dict /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x.ground_truth.tsv --output /cluster/dataset/nexus/benchmarking/celltyping/data/summary_tables/ --sample_name Adul
```

Output directory: /cluster/dataset/nexus/benchmarking/celltyping/data/summary_tables/

##########################################

## (e) Summary of paths to relevant Adult PBMC files

- Barcodes: /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x.scp.barcodes.tsv
- Metadata: /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x_annotated_v1.scp.metadata.tx
- Ground-truth: /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Adult-Immune-10x/2020-Mar-Census-Adult-Immune-10x.ground_truth.tsv

- Processing run directory: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_scrun/
- Processing run (subset) directory: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_subset_scrun/

- Data hdf5 (raw counts): /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_scrun/analysis/rawCounts/2020-Mar-Census-Adult-Immune-10x.h5
- Subset hdf5 (raw counts): /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/2020-Mar-Census-Adult-Immune-10x.scp_subset.h5

- Processed RDS: /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Adult-Immune-10x_run/analysis/prep_celltyping/2020-Mar-Census-Adult-Immune-10x.genes_cells_filtered.corrected.RDS #on Euler
- Processed RDS (subset): /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020_09_18-Zheng_subset_scrun/analysis/prep_celltyping/Zheng_subset.genes_cells_filtered.corrected.RDS

- Annotated RDS: [N/A] - has not been created.
- Annotated RDS (subset): /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/2020-Mar-Census-Adult-Immune-10x_subset_annotated.RDS

- Proportions table: /cluster/work/nexus/benchmarking/celltyping/data/summary_tables/adult/Adult_true_labels_table.tx

----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------

# [3] 2020-Mar-Census-Newborn-Blood-10x datase
Newborn PBMC dataset (preprocessed) retrieved from the Broad Institute.
See /cluster/work/nexus/benchmarking/celltyping/data/broad/readme.txt for more details.

Since the data is already preprocessed (cellranger output), we can proceed directly to the processing step.

----------------------------------------------------------------------------------------------------------------
## (a) Processing - running scRNAseq pipeline on dataset

WARNING: full dataset requires a special implementation of the processing pipeline due to the size of the dataset.
As such, only the subset will be used for annotation.

Commands are still experimental and all files can be found on Euler.

Snakemake command used (on Euler):
```
export HDF5_USE_FILE_LOCKING=FALSE; /cluster/project/nexus/utilities/sharedPrograms/snakemake/snakemake_4.8.0/bin/snakemake --notemp --latency-wait 240 -s /cluster/project/nexus/benchmarking/celltyping/test_pilot/2020-Mar-Census-Newborn-Blood-10x_run/snake_analysis_files/snake_gex_master.snake --configfile /cluster/project/nexus/benchmarking/celltyping/test_pilot/2020-Mar-Census-Newborn-Blood-10x_run/snake_analysis_files/config.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R rusage[mem={params.mem},scratch={params.scratch}] -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k;
```

Output location: /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_run

##########################################

## (b) Generating subset

The first step of the scRNAseq processing pipeline is generating an hdf5 file from the cellranger count output containing the raw count data.
This is used to generate an hdf5 file containing a subset of the original data. Around 10% of cells in the original dataset are chosen randomly, respecting and maintaing the proportion of cell types in the dataset.

Subsetting was performed with an R script as follows (on Euler):
`Rscript /cluster/project/nexus/marcus/git/celltype_benchmark/scripts/subset_data_proportional.R --input_data /cluster/project/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_run/analysis/rawCounts/2020-Mar-Newborn-Blood-10x.h5 --true_labels /cluster/project/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-Blood-10x.ground_truth.tsv --output_subset /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/2020-Mar-Census-Newborn-Blood-10x.scp_subset.h5`

Output location: /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/

----------------------------------------------------------------------------------------------------------------
## i. Processing - running scRNAseq pipeline on subset

Snakemake command used (on Euler):
```
export HDF5_USE_FILE_LOCKING=FALSE; /cluster/project/nexus/utilities/sharedPrograms/snakemake/snakemake_4.8.0/bin/snakemake --notemp --latency-wait 240 -s /cluster/project/nexus/benchmarking/celltyping/test_pilot/2020-Mar-Census-Adult-Immune-10x_run/snake_analysis_files/snake_gex_master.snake --configfile /cluster/project/nexus/benchmarking/celltyping/test_pilot/newborn_subset/snake_analysis_files/config.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R rusage[mem={params.mem},scratch={params.scratch}] -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k;
```

Output location: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_subset_scrun

Prior to the run, a symbolic link to the subset hdf5 file is created in the pipeline run directory:
`ln -s /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/2020-Mar-Census-Newborn-Blood-10x.scp_subset.h5 /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_subset_scrun/analysis/rawCounts/2020-Mar-Census-Newborn-Blood-10x_subset.scp.h5`

##########################################

## (c) Annotating data with cell type ground-truth labels

For the benchmark, the data must first be labeled with the true cell types per barcode.
To do this, a file (hereafter referred to as a ground-truth dictionary) mapping the labels to each barcode must be generated.
Note that this only needs to be done once per (full) data set.
Please see the documentation on the celltyping_benchmark Git repo for examples and required formatting of the input files needed.

Command used to generate ground-truth dictionary for the Adult PBMC data:
`python /cluster/work/nexus/benchmarking/celltyping/git/scripts/generate_ground-truth_list.py --original_labels /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-Blood-10x_annotated_v1.scp.metadata.txt --label_dictionary /cluster/work/nexus/benchmarking/celltyping/git/doc/matched_types_dictionary.tsv --output /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-Blood-10x.ground_truth.tsv --subtype_dictionary /cluster/work/nexus/benchmarking/celltyping/data/major_subtype_mapping.txt`

----------------------------------------------------------------------------------------------------------------

Once this is done, the processed RDS file from the scRNAseq processing pipeline's 'prep_celltyping' step is ready to be annotated.

Command used to annotate processed RDS with ground-truth labels (Newborn PBMC subset):
```
Rscript /cluster/work/nexus/benchmarking/celltyping/git/scripts/add_sce_groundtruth.R --sce_in cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_subset_scrun/analysis/prep_celltyping/2020-Mar-Census-Newborn-Blood-10x_subset.scp.genes_cells_filtered.corrected.RDS --sce_out /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/2020-Mar-Census-Newborn-Blood-10x_subset_annotated.RDS --labels_in /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-10x_annotated_v1.scp.metadata.tx
```

Output path: /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/

##########################################

## (d) Table of cell type proportions in dataset (full & subset) 

A script was used to generated tables summarizing each dataset's cell type counts and proportions.
The prefixes of the column labels map as follows:
- "raw_full": barcodes file (cellranger output)
- "raw_subset": HDF5 file of subset data (file generated by subsetting script)
- "processed_full": RDS file of SCE object (output of scRNAseq pipeline: 'prep_celltyping')
- "processed_subset": RDS file of subsetted SCE object (output of scRNAseq pipeline: 'prep_celltyping')
- "metadata": ground-truth dictionary
("processed_full" omitted due to size of file/object.)

Command used to generated Newborn PBMC summary table (on Leomed):
```
Rscript /cluster/work/nexus/benchmarking/celltyping/data/summary_tables/write_celltype_proportion_tables.R --data_raw /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_scrun/analysis/rawCounts/2020-Mar-Census-Newborn-Blood-10x.h5 --data_subset /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/2020-Mar-Census-Newborn-Blood-10x_subset_annotated.RDS --data_subset_raw /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/2020-Mar-Census-Newborn-Blood-10x.scp_subset.h5 --data_dict /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-Blood-10x.ground_truth.tsv --output /cluster/dataset/nexus/benchmarking/celltyping/data/summary_tables/ --sample_name Newborn
```

Output directory: /cluster/dataset/nexus/benchmarking/celltyping/data/summary_tables/

##########################################

## (e) Summary of paths to relevant Newborn PBMC files (on Leomed)

- Barcodes: /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-Blood-10x.scp.barcodes.tsv
- Metadata: /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-10x_annotated_v1.scp.metadata.tx
- Ground-truth: /cluster/work/nexus/benchmarking/celltyping/data/broad/2020-Mar-Census-Newborn-Blood-10x/2020-Mar-Census-Newborn-Blood-10x.ground_truth.tsv

- Processing run directory: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_scrun/
- Processing run (subset) directory: /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_subset_scrun/

- Data hdf5 (raw counts): /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_scrun/analysis/rawCounts/2020-Mar-Census-Newborn-Blood-10x.h5
- Subset hdf5 (raw counts): /cluster/work/nexus/benchmarking/celltyping/data/subset/raw_counts_hdf5/2020-Mar-Census-Newborn-Blood-10x.scp_subset.h5

- Processed RDS: [N/A] - has not been created.
- Processed RDS (subset): /cluster/work/nexus/benchmarking/celltyping/preprocessing_runs/2020-Mar-Census-Newborn-Blood-10x_subset_scrun/analysis/prep_celltyping/2020-Mar-Census-Newborn-Blood-10x_subset.scp.genes_cells_filtered.corrected.RDS

- Annotated RDS: [N/A] - has not been created.
- Annotated RDS (subset): /cluster/work/nexus/benchmarking/celltyping/data/subset/extended_sce/2020-Mar-Census-Newborn-Blood-10x_subset_annotated.RDS

- Proportions table: /cluster/work/nexus/benchmarking/celltyping/data/summary_tables/newborn/Newborn_true_labels_table.tx


