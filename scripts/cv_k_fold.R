#!/usr/bin/env Rscript

################################################################################
## creating the k-fold crossvalidation data
################################################################################
#library
library(optparse)
library(tidyverse)
library(caret)

#Data Path

# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_basic.RDS)."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--kfold", type = "integer", help = " k-fold cross-validation: The original sample is partitioned into k equal sized subsamples."),
  make_option("--sample_name", type = "character", help = "Name of sample to be appended to all output files.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

################################################################################
## main code starts here
################################################################################
## load input data
sce_data <- readRDS(opt$SCE)
lab_maj <- sce_data@metadata$ground_truth_major
lab_data <- sce_data@metadata$ground_truth_minor
lab_data[lab_data=="none"] <- lab_maj[lab_data=="none"]
barcode <- colnames(sce_data)

# cross validation
##################
folds = createFolds(lab_data, k = opt$kfold,list = T)
fold_names <- str_split(names(folds),"d",simplify = T)
fold_names <- apply(fold_names,1,function(x){paste(x[1],"d_",x[2],sep = "")})

allIndex <- 1:length(lab_data)
for(i in 1:length(folds)){
  test_idx <- folds[[i]]
  train_idx <- allIndex[-folds[[i]]]
  test_data <- barcode[test_idx]
  train_data <- barcode[train_idx]
  saveRDS(list(test_data=test_data,train_data=train_data),
            paste(opt$outputDirec,opt$sample_name,".index",fold_names[i],".RDS",sep=""))
}
saveRDS(list(test_data=numeric(),train_data=barcode[allIndex]),paste(opt$outputDirec,opt$sample_name,".index_All.RDS",sep=""))


