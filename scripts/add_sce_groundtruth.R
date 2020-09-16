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
    sce@metadata$true_labels <- runif(ncol(sce))
    sce@metadata$ground_truth_major <- runif(ncol(sce))
    sce@metadata$ground_truth_minor <- runif(ncol(sce))

    labels <- read.table(labels_path, header = FALSE, sep ='\t', skip = 1)
    true_major_label <- labels[,2]
    print(head(true_major_label))
    matched_major_label <- labels[,3]
    matched_minor_label <- labels[,4]
    names(true_major_label) <- labels[,2]
    names(matched_major_label) <- labels[,3]
    names(matched_minor_label) <- labels[,4]
    
    #match_index <- match(sce$barcodes, substr(names(matched_major_label),1,16))
    match_index <- match(sce$barcodes, labels[,1])
    if(NA %in% match_index){
	unmatched_index <- which(is.na(match_index))
        print(unmatched_index)
	missing_index <- match(sce$barcodes[unmatched_index], substr(labels[,1],1,16))
	print(sce$barcodes[unmatched_index]);print(substr(labels[,1][unmatched_index],1,16));
	match_index[unmatched_index] <- missing_index
    }
    if(NA %in% match_index){
        print(paste("Warning:",length(which(is.na(match_index))),"cell barcodes not matched to cell types."))
    }
    sce@metadata$ground_truth_major <- as.character(matched_major_label[match_index])
    sce@metadata$ground_truth_minor <- as.character(matched_minor_label[match_index])
    sce@metadata$true_labels <- as.character(true_major_label[match_index])
    
    saveRDS(sce, output_path)
    }

add_sce_groundtruth(opt$sce_in, opt$labels_in, opt$sce_out)
