#### Garnett - training script ####

# Training component of Garnett cell typing method.
# Read in processed SCE & marker genes file to train a classifier on provided counts data and markers.

library(garnett)
library(org.Hs.eg.db)

option_list = list(
  make_option("--cds", type = "character", help = "Path to RDS file with CDS data object."),
  make_option("--barcodes_index", type = "character", help = "Path to RDS file w/ train & test data barcodes."),
  make_option("--output_dir", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sample_name", type = "character", help = "Sample identifier. Attached to each output name."),
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cds <- readRDS(opt$cds)
barcodes_selected <- readRDS(opt$barcodes_index)
cds <- estimate_size_factors(cds)

######### Garnett run ##########    
print("Starting training...")
start_train <- Sys.time()
garnett_classifier <- train_cell_classifier(cds = cds[,barcodes_selected$train_data], 
                                           marker_file = opt$marker_file,
                                           db=org.Hs.eg.db,
                                           cds_gene_id_type = "SYMBOL",
                                           num_unknown = 50,
                                           marker_file_gene_id_type = "SYMBOL")
end_train <- Sys.time()
print("Training complete.")
train_time <- as.numeric(end_train - start_train)
write.csv(train_time,paste0(opt$output_dir, opt$sample_name,'/garnett_training_time.csv'),row.names = FALSE)

print("Testing classifier; generating predicted labels...")
start_test <- Sys.time()
garnett_test <- classify_cells(cds[,barcodes_selected$test_data],
                                  garnett_classifier,
                                  db = org.Hs.eg.db,
                                  cluster_extend = TRUE,
                                  cds_gene_id_type = "SYMBOL")

end_test <- Sys.time()
test_time <- as.numeric(end_test - start_test)
print("Prediction complete.")
print("Saving predicted labels...")
pred_labels <- list(pData(cds)$cluster_ext_type)
write.csv(pred_labels,paste0(opt$output_dir, opt$sample_name, '/garnett_pred_label.csv'),row.names = FALSE)
write.csv(test_time,paste0(opt$output_dir, opt$sample_name, '/garnett_CV_test_time.csv'),row.names = FALSE)
print("Labels saved.")

print("Saving classifier model...")
saveRDS(garnett_classifier, paste0(opt$output_dir, opt$sample_name, ".garnett_model.RDS"))
print("Model saved.")
