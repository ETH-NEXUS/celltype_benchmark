#### Garnett - testing script ####

# Testing component of Garnett cell typing method.
# Read in pre-trained classifier & a processed SCE file containing counts data to classify cells by cell types identified in training step.

library(garnett)
library(org.Hs.eg.db)

option_list = list(
  make_option("--sce", type = "character", help = "Path to RDS file with sce object stored inside."),
  make_option("--output_dir", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sample_name", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--classifier", type="character", help = "Path to pre-trained Garnett classifier generated in training step.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

sce <- readRDS(opt$sce)
sce <- estimate_size_factors(sce)

garnett_classifier <- readRDS(opt$classifier)
    
print("Starting testing...")
start_test <- Sys.time()
  
garnett_test <- classify_cells(sce, 
                                  garnett_classifier, 
                                  db = org.Hs.eg.db, 
                                  cluster_extend = TRUE,
                                  cds_gene_id_type = "SYMBOL")

end_test <- Sys.time()
test_time <- as.numeric(end_test - start_test)
  
#true_labels <- list(sce$true_label_major)
pred_labels <- list(pData(sce)$cluster_ext_type)

print("Writing test results...")
#write.csv(true_labels,paste0(opt$output_dir,'/Garnett_CV_true.csv'),row.names = FALSE)
write.csv(pred_labels,paste0(opt$output_dir,'/Garnett_CV_pred.csv'),row.names = FALSE)
write.csv(test_time,paste0(opt$output_dir,'/Garnett_CV_test_time.csv'),row.names = FALSE)
print("Done!")
