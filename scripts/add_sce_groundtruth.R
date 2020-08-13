library(optparse)
library(SingleCellExperiment)

option_list = list(
  make_option("--sce_in", type = "character", help = "Input path to processed SCE object (.RDS format)."),
  make_option("--labels_in", type = "character", help = "Path to ground-truth label file. Required structure: column 1 = barcodes (full length), column 2 = major cell type, column 3 = minor cell type."),
  make_option("--sce_out", type = "character", help = "Output path of annotated SCE object (.RDS format).")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

add_sce_groundtruth <- function(rds_path, labels_path, output_path){
    sce <- readRDS(rds_path)
    sce@metadata$ground_truth_major <- runif(ncol(sce))
    sce@metadata$ground_truth_minor <- runif(ncol(sce))

    labels <- read.table(labels_path, header = TRUE, sep ='\t')
    major_label <- labels$major_label
    minor_label <- labels$minor_label
    names(major_label) <- labels$barcode
    names(minor_label) <- labels$barcode
    
    match_index <- match(sce$barcodes, substr(names(major_label),1,16))

    sce@metadata$ground_truth_major <- as.character(major_label[match_index])
    sce@metadata$ground_truth_minor <- as.character(minor_label[match_index])
    
    saveRDS(sce, output_path)
    }

add_sce_groundtruth(opt$sce_in, opt$labels_in, opt$sce_out)
