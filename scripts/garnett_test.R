#### Garnett - testing script ####

# Testing component of Garnett cell typing method.
# Read in pre-trained classifier & a processed SCE file containing counts data to classify cells by cell types identified in training step.

library(optparse)
library(garnett)
library(org.Hs.eg.db)

option_list = list(
  make_option("--cds", type = "character", help = "Path to RDS file with cds object stored inside."),
  make_option("--outputFile", type = "character", help = "Path to the output file."),
#  make_option("--sample_name", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--classifier", type="character", help = "Path to pre-trained Garnett classifier generated in training step.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cds <- readRDS(opt$cds)
cds <- estimate_size_factors(cds)

garnett_classifier <- readRDS(opt$classifier)
    
print(paste("Sample:", opt$cds))
print(paste("Model:",  opt$classifier))

print("Starting testing...")
start_test <- Sys.time()
  
cds <- classify_cells(cds, 
                      garnett_classifier, 
                      db = org.Hs.eg.db, 
                      cluster_extend = TRUE,
                      cds_gene_id_type = "ENSEMBL",
		      verbose = TRUE)

end_test <- Sys.time()
test_time <- as.numeric(end_test - start_test)
  
#true_labels <- list(cds$true_label_major)
pred_labels <- list(pData(cds)$cluster_ext_type)

print("Writing predicted labels...")
#write.csv(true_labels,paste0(opt$output_dir,'/Garnett_CV_true.csv'),row.names = FALSE)
#write.csv(pred_labels,opt$outputFile,row.names = FALSE)
library(tidyverse)
write_csv(data.frame(celltype_final=pred_labels[[1]]),opt$outputFile)
#write.csv(pred_labels,paste0(opt$output_dir, opt$sample_name, '.garnett_CV_pred.csv'),row.names = FALSE)
#write.csv(test_time,paste0(opt$output_dir, opt$sample_name, '.garnett_CV_test_time.csv'),row.names = FALSE)
print("Done!")
