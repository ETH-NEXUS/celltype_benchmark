Description of data set processing for the scROSHI publication

Note: 
/PATH_TO_GIT/ refers to the path with the celltype benchmark repository
/PATH_TO_PROJECT/ refers (for the purpose of this documentation) to the directory with all scROSHI relevant data (no public access) 

# Example data - processing workflows

There are 3 data sets provided for an example run of the cell typing benchmark.
 
- Zheng (sorted & merged)
- (2020 Broad) Adult PBMC
- (2020 Broad) Newborn PBMC

Their pre-processing, processing, and annotation steps are detailed below.

######################
Pre-processing and processing necessary for each respective datase
#####################

## [1] Zheng data set
Consists of 10 FACS sorted immune cell populations which are then merged into one dataset.
Data sets:
1. CD19+ B cells
2. CD14+ Monocytes
3. CD34+ Cells
4. CD4+ Helper T cells
5. CD56+ Natural Killer Cells
6. CD8+ Cytotoxic T cells
7. CD4+/CD45RO+ Memory T cells
8. CD8+/CD45RA+ Naive Cytotoxic T cells
9. CD4+/CD45RA+/CD25- Naive T cells
10. CD4+/CD25+ Regulatory T cells

FASTQ files for each dataset were retrieved from the '10X Genomics Single Cell Gene Expression Datasets' page.
(https://support.10xgenomics.com/single-cell-gene-expression/datasets - Single Cell 3' Paper: Zheng et al. 2017)

----------------------------------------------------------------------------------------------------------------
### (a) Pre-processing - cellranger count

Each dataset was processed separately with the cellranger count pipeline (cellranger version 3.1.0), mapped against the GRCh38 reference.

Example call for the CD19+ B cells:

```
cellranger count --id=b_cells --transcriptome=/PATH_TO_PROJECT/refdata-cellranger-GRCh38-3.0.0 --localcores=24 --fastqs=/PATH_TO_PROJECT/b_cells_fastqs --nosecondary
```

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
b_cells,molecule_info.h5,b_cells
```

Command to merge: 
`cellranger aggr --id=zheng_merged --csv=/PATH_TO_PROJECT/aggregate_libraries.csv --normalize=none`

----------------------------------------------------------------------------------------------------------------
### (b) Processing of cellranger aggr output into an h5 file (using the scAmpi pipeline (Bertolini et al., 2021))
(Note: only the basic step of creating the h5 file from the raw cellranger output was of interest)

Basic snakemake command (details provided in scAmpi documentation):
```
snakemake -s scAmpi_gex_master.snake --configfile config_scAmpi.json -j 100 -p -k
```

##########################################

### (c) Generating subset and creating an SCE
The pooled data sets used in the benchmark study contain several thousands of cells, which is a significantly larger number compared to a typical 10xGenomics scRNA-seq data set. To avoid resource issued while preserving a data set size comparable with typical experiments, we downsampled the original data sets to 10% of their size. 
Around 10% of cells in the original dataset are chosen randomly, respecting and maintaing the proportion of cell types in the dataset.

Subsetting was performed with an R script as follows (on the h5 file generated by the scAmpi pipeline):
`Rscript /PATH_TO_GIT/scripts/subset_data_proportional.R --input_data /PATH_TO_PROJECT/Zheng_merged.h5 --true_labels /PATH_TO_PROJECT/zheng_ground_truth.tsv --output_subset /PATH_TO_PROJECT/Zheng_subset.h5`

The subset h5 file was processed into an SCE object using the scAmpi pipeline (Bertolini et al., 2021)
(Note: only basic contamination removal and count normalization were of interest, clustering and cell typing can be ignored).
----------------------------------------------------------------------------------------------------------------

##########################################

### (d) Annotating data with cell type ground-truth labels

For the benchmark, the data must first be labeled with the true cell types per barcode.
To do this, first a file (hereafter referred to as a ground-truth dictionary) mapping the labels to each barcode must be generated.

Command used to generate ground-truth dictionary for the Zheng data:
`python /PATH/TO/GIT/scripts/generate_ground-truth_list.py --original_labels /PATH_TO_PROJECT/zheng_original.tsv --label_dictionary /PATH/TO/GIT/doc/example_matched_types_dictionary.tsv --output /PATH_TO_PROJECT/zheng_ground_truth.tsv --subtype_dictionary /PATH/TO/GIT/doc/doc_scROSHI_paper/major_subtype_mapping.txt`

Now the RDS file generated by scAmpi (containing the SCE object with count data) can be annotated  with the ground-truth labels (Zheng data):
```
Rscript /PATH/TO/GIT/scripts/add_sce_groundtruth.R --sce_in /PATH_TO_PROJECT/Zheng_subset.RDS --sce_out /PATH_TO_PROJECT/Zheng_merged_annotated.RDS --labels_in /PATH_TO_PROJECT/zheng_ground_truth.tsv
```
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------

## [2] 2020-Mar-Census-Adult-Immune-10x datase
Adult PBMC dataset (preprocessed) retrieved from the Broad Institute.
See /cluster/work/nexus/benchmarking/celltyping/data/broad/readme.txt for more details.

Since the data is already preprocessed (cellranger output), we can proceed directly to the processing step.

----------------------------------------------------------------------------------------------------------------
### (a) Processing of cellranger output into an h5 file (using the scAmpi pipeline (Bertolini et al., 2021))
(Note: only the basic step of creating the h5 file from the raw cellranger output was of interest)

Basic snakemake command (details provided in scAmpi documentation):
```
snakemake -s scAmpi_gex_master.snake --configfile config_scAmpi.json -j 100 -p -k
```

##########################################

### (b) Generating subset and creating an SCE
The pooled data sets used in the benchmark study contain several thousands of cells, which is a significantly larger number compared to a typical 10xGenomics scRNA-seq data set. To avoid resource issued while preserving a data set size comparable with typical experiments, we downsampled the original data sets to 10% of their size. 
Around 10% of cells in the original dataset are chosen randomly, respecting and maintaing the proportion of cell types in the dataset.

Subsetting was performed with an R script as follows (on the h5 file generated by the scAmpi pipeline):
`Rscript /PATH_TO_GIT/scripts/subset_data_proportional.R --input_data /PATH_TO_PROJECT/2020-Mar-Census-Adult-Immune-10x.h5 --true_labels /PATH_TO_PROJECT/2020-Mar-Census-Adult-Immune-10x.ground_truth.tsv --output_subset /PATH_TO_PROJECT/2020-Mar-Census-Adult-Immune-10x.scp_subset.h5`

The subset h5 file was processed into an SCE object using the scAmpi pipeline (Bertolini et al., 2021)
(Note: only basic contamination removal and count normalization were of interest, clustering and cell typing can be ignored).

##########################################

### (c) Annotating data with cell type ground-truth labels

For the benchmark, the data must first be labeled with the true cell types per barcode.
To do this, first a file (hereafter referred to as a ground-truth dictionary) mapping the labels to each barcode must be generated.

Command used to generate ground-truth dictionary for the Zheng data:
`python /PATH/TO/GIT/scripts/generate_ground-truth_list.py --original_labels /PATH_TO_PROJECT/2020-Mar-Census-Adult-Immune-10x_annotated_v1.scp.metadata.txt --label_dictionary /PATH/TO/GIT/doc/example_matched_types_dictionary.tsv --output /PATH_TO_PROJECT/2020-Mar-Census-Adult-Immune-10x.ground_truth.tsv --subtype_dictionary /PATH/TO/GIT/doc/doc_scROSHI_paper/major_subtype_mapping.txt`

Now the RDS file generated by scAmpi (containing the SCE object with count data) can be annotated  with the ground-truth labels (Zheng data):
```
Rscript /PATH/TO/GIT/scripts/add_sce_groundtruth.R --sce_in /PATH_TO_PROJECT/2020-Mar-Census-Adult-Immune-10x_subset.RDS --sce_out /PATH_TO_PROJECT/2020-Mar-Census-Adult-Immune-10x_subset_annotated.RDS --labels_in /PATH_TO_PROJECT/2020-Mar-Census-Adult-Immune-10x.ground_truth.tsv
```

----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------

# [3] 2020-Mar-Census-Newborn-Blood-10x datase
Newborn PBMC dataset (preprocessed) retrieved from the Broad Institute.
See /cluster/work/nexus/benchmarking/celltyping/data/broad/readme.txt for more details.

Since the data is already preprocessed (cellranger output), we can proceed directly to the processing step.

----------------------------------------------------------------------------------------------------------------
### (a) Processing of cellranger output into an h5 file (using the scAmpi pipeline (Bertolini et al., 2021))
(Note: only the basic step of creating the h5 file from the raw cellranger output was of interest)

Basic snakemake command (details provided in scAmpi documentation):
```
snakemake -s scAmpi_gex_master.snake --configfile config_scAmpi.json -j 100 -p -k
```

##########################################

### (b) Generating subset and creating an SCE
The pooled data sets used in the benchmark study contain several thousands of cells, which is a significantly larger number compared to a typical 10xGenomics scRNA-seq data set. To avoid resource issued while preserving a data set size comparable with typical experiments, we downsampled the original data sets to 10% of their size. 
Around 10% of cells in the original dataset are chosen randomly, respecting and maintaing the proportion of cell types in the dataset.

Subsetting was performed with an R script as follows (on the h5 file generated by the scAmpi pipeline):
`Rscript /PATH_TO_GIT/scripts/subset_data_proportional.R --input_data /PATH_TO_PROJECT/2020-Mar-Newborn-Blood-10x.h5 --true_labels /PATH_TO_PROJECT/2020-Mar-Census-Newborn-Blood-10x.ground_truth.tsv --output_subset /PATH_TO_PROJECT/2020-Mar-Census-Newborn-Blood-10x.scp_subset.h5`

The subset h5 file was processed into an SCE object using the scAmpi pipeline (Bertolini et al., 2021)
(Note: only basic contamination removal and count normalization were of interest, clustering and cell typing can be ignored).

##########################################

### (c) Annotating data with cell type ground-truth labels

For the benchmark, the data must first be labeled with the true cell types per barcode.
To do this, first a file (hereafter referred to as a ground-truth dictionary) mapping the labels to each barcode must be generated.

Command used to generate ground-truth dictionary for the Zheng data:
`python /PATH/TO/GIT/scripts/generate_ground-truth_list.py --original_labels /PATH_TO_PROJECT/2020-Mar-Census-Newborn-Blood-10x_annotated_v1.scp.metadata.txt --label_dictionary /PATH/TO/GIT/doc/example_matched_types_dictionary.tsv --output /PATH_TO_PROJECT/2020-Mar-Census-Newborn-Blood-10x.ground_truth.tsv --subtype_dictionary /PATH/TO/GIT/doc/doc_scROSHI_paper/major_subtype_mapping.txt`

Now the RDS file generated by scAmpi (containing the SCE object with count data) can be annotated  with the ground-truth labels (Zheng data):
```
Rscript /PATH/TO/GIT/scripts/add_sce_groundtruth.R --sce_in /PATH_TO_PROJECT/2020-Mar-Census-Newborn-Blood-10x_subset.RDS --sce_out /PATH_TO_PROJECT/2020-Mar-Census-Newborn-Blood-10x_subset_annotated.RDS --labels_in /PATH_TO_PROJECT/2020-Mar-Census-Newborn-Blood-10x.ground_truth.tsv
```
