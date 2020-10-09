#### Garnett - training script ####

# Training component of Garnett cell typing method.
# Read in processed SCE & marker genes file to train a classifier on provided counts data and markers.

library(optparse)
library(garnett)
library(org.Hs.eg.db)

option_list = list(
  make_option("--cds", type = "character", help = "Path to RDS file with CDS data object."),
  make_option("--barcodes_index", type = "character", help = "Path to RDS file w/ train & test data barcodes."),
  make_option("--output_dir", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sample_name", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--marker_file", type = "character", help = "Path to Garnett-formatted marker gene list.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
set.seed(260)
cores_to_use <- parallel::detectCores() - 1

cds <- readRDS(opt$cds)
barcodes_selected <- readRDS(opt$barcodes_index)
cds <- estimate_size_factors(cds)

# Check gene ID type of genes in marker file & CDS
example_gene_in_marker_file <- strsplit(readLines(opt$marker_file, n = 2)[2],", ")[[1]][2]
marker_gene_format <- ifelse(substr(example_gene_in_marker_file, start = 1, stop = 4) == "ENSG", "ENSEMBL", "SYMBOL")

example_gene_in_cds <- rownames(cds)[1]
cds_gene_format <- ifelse(substr(example_gene_in_cds, start = 1, stop = 4) == "ENSG", "ENSEMBL", "SYMBOL")

print(paste("Gene name format in marker file:", marker_gene_format))
print(paste("Gene name format in CDS data:", cds_gene_format))

######### Garnett run ##########    
print("Starting training...")
start_train <- Sys.time()
garnett_classifier <- train_cell_classifier(cds = cds[,barcodes_selected$train_data], 
                                           marker_file = opt$marker_file,
                                           db=org.Hs.eg.db,
                                           cds_gene_id_type = cds_gene_format,
					   max_training_samples = dim(cds)[2],
					   min_observations = 20,
                                           num_unknown = 500,
					   propogate_markers = TRUE,
                                           marker_file_gene_id_type = marker_gene_format,
					   cores = cores_to_use)
end_train <- Sys.time()
print("Training complete.")
train_time <- as.numeric(end_train - start_train)
#write.csv(train_time,paste0(opt$output_dir, opt$sample_name,'.garnett_training_time.csv'),row.names = FALSE)

if(length(barcodes_selected$test_data) > 0){
    print("Testing classifier; generating predicted labels...")
    start_test <- Sys.time()
    garnett_test <- classify_cells(cds[,barcodes_selected$test_data],
                                      garnett_classifier,
                                      db = org.Hs.eg.db,
                                      cluster_extend = TRUE,
                                      cds_gene_id_type = cds_gene_format)

    end_test <- Sys.time()
    test_time <- as.numeric(end_test - start_test)
    print("Prediction complete.")
    print("Saving predicted labels...")
    print(paste0(opt$output_dir, opt$sample_name, '.garnett_predicted_labels.csv'))
    pred_labels <- list(pData(garnett_test)$cluster_ext_type)
    names(pred_labels) <- "x"
    write.csv(pred_labels,paste0(opt$output_dir, opt$sample_name, '.garnett_predicted_labels.csv'),row.names = FALSE,quote=FALSE)
    #write.csv(test_time,paste0(opt$output_dir, opt$sample_name, '.garnett_CV_test_time.csv'),row.names = FALSE)
    print("Labels saved.")
}

print("Saving classifier model...")
print(paste0(opt$output_dir, opt$sample_name, ".garnett_model.RDS"))
saveRDS(garnett_classifier, paste0(opt$output_dir, opt$sample_name, ".garnett_model.RDS"))
print("Model saved.")

debug_markers  <- train_cell_classifier(cds = cds[,barcodes_selected$train_data],
                                           marker_file = opt$marker_file,
                                           db=org.Hs.eg.db,
                                           cds_gene_id_type = cds_gene_format,
                                           max_training_samples = dim(cds)[2],
                                           min_observations = 20,
                                           num_unknown = 50,
                                           propogate_markers = TRUE,
                                           marker_file_gene_id_type = marker_gene_format,
                                           cores = 1, return_initial_assign = TRUE)

print(debug_markers)
