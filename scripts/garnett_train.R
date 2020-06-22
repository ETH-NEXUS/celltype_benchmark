#### Garnett - training script ####

# Training component of Garnett cell typing method.
# Read in processed SCE & marker genes file to train a classifier on provided counts data and markers.

library(garnett)
library(org.Hs.eg.db)

option_list = list(
  make_option("--sce", type = "character", help = "Path to RDS file with sce object stored inside."),
  make_option("--output_dir", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sample_name", type = "character", help = "Sample identifier. Attached to each output name."),
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

sce <- readRDS(opt$sce)
sce <- estimate_size_factors(sce)
    
print("Starting training...")
start_train <- Sys.time()
garnett_classifier <- train_cell_classifier(cds = sce, 
                                           marker_file = opt$marker_file,
                                           db=org.Hs.eg.db,
                                           cds_gene_id_type = "SYMBOL",
                                           num_unknown = 50,
                                           marker_file_gene_id_type = "SYMBOL")
end_train <- Sys.time()
train_time <- as.numeric(end_train - start_train)

write.csv(train_time,paste0(opt$output_dir,'/garnett_training_time.csv'),row.names = FALSE)
saveRDS(garnett_classifier, paste0(opt$output_dir, "garnett_classifier.RDS"))
