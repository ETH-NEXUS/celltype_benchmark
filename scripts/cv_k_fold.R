################################################################################
## creating the k-fold crossvalidation data
################################################################################

#library
library(optparse)
library(tidyverse)
library(caret)

#Data Path
#opt = list(
#  outputDirec = "/Users/bolars/Documents/celltyping/benchmark_scripts/",
#  SCE = "/Users/bolars/Documents/celltyping/benchmark_scripts/Zheng_sorted_merged.genes_cells_filtered.corrected.ground-truth.RDS"
#)

# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_basic.RDS)."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Parameters
Kfold <- 10

################################################################################
## main code starts here
################################################################################
## load input data
sce_data <- readRDS(opt$SCE)
lab_data <- colData(sce_data)$true_label
barcode <- colnames(sce_data)

# cross validation
##################
folds = createFolds(lab_data, k = Kfold,list = T)
allIndex <- 1:length(lab_data)
for(i in 1:length(folds)){
  test_idx <- folds[[i]]
  train_idx <- allIndex[-folds[[i]]]
  test_data <- barcode[test_idx]
  train_data <- barcode[train_idx]
  saveRDS(list(test_data=test_data,train_data=train_data),
            paste(opt$outputDirec,"index",names(folds)[i],".RDS",sep=""))
}
saveRDS(list(test_data=numeric(),train_data=barcode[allIndex]),paste(opt$outputDirec,"indexAll.RDS",sep=""))

