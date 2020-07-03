# Script for converting SCE object (saved as RDS) into a Monocle3-compliant Cell_Data_Set (CDS) object as input for the Garnett cell typer.
# Requires conda environment installed via /doc/conda_envs/garnett.yaml

library(monocle3)
library(org.Hs.eg.db)

option_list = list(
  make_option("--sce_in", type = "character", help = "Path to RDS file containing SCE data object."),
  make_option("--cds_out", type = "character", help = "Path to output file of final (converted) CDS object."),
  make_option("--sample_name", type = "character", help = "Sample name of input SCE data.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Read in SCE object & list of barcodes to subset it on
data <- readRDS(opt$sce_in)

# Preserve data/annotations that are not automatically included in conversion
umap_cl <- data@colData$umap_cl
umap_hvg <- reducedDim(data) #cols: barcode, X, Y?
ground_truth_labels <- c(data@metadata$ground_truth_major, data@metadata$ground_truth_minor)
gene_ids <- rowData(data)$gene_ids
normcounts <- data@assays$data$normcounts

# Create a CDS object using counts data
data <- monocle3::new_cell_data_set(expression_data=data@assays$data$counts)

# Add data/annotations manually to converted CDS object
fData(data)$gene_id <- gene_ids
fData(data)$umap_cl <- umap_cl
fData(data)$ground_truth_major <- ground_truth_labels[1]
fData(data)$ground_truth_minor <- ground_truth_labels[2]
reducedDims(data)$UMAP <- umap_hvg

# Add shortened gene names to CDS
features <- read.table("/cluster/project/nexus/benchmarking/celltyping/test_pilot/cellranger_run/zheng_sorted_merged/outs/filtered_feature_bc_matrix/features.tsv.gz") #TO-DO: add universal features files in a /data/ dir & point to it using config file.
features <- features[1:2]
names(features) <- c("gene_id","gene_short_name")
match_index <- match(row.names(data), features$gene_id)
fData(data)$gene_short_name <- features$gene_short_name[match_index]

#Save newly converted CDS file(s): train set & test set for Garnett classifier training.
saveRDS(data, paste0(opt$output_dir, opt$sample_name, ".monocle3_cds.RDS"))
