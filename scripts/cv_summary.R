################################################################################
## Summary of the k-fold crossvalidation
################################################################################
#library
library(optparse)
library(tidyverse)
library(e1071)
library(cowplot)
library(caret)
library(Hmisc)

#Data Path
#opt = list(
#  outputDirec = "/Users/bolars/Documents/celltyping/benchmark_scripts/",
#  CVfiles = "/Users/bolars/Documents/celltyping/benchmark_scripts/",
#  SCE = "/Users/bolars/Documents/celltyping/benchmark_scripts/Zheng_sorted_merged.genes_cells_filtered.corrected.ground-truth.RDS",
#  sampleName = "CV_summary",
#  method = c("RF","SVM","fcNN")
#)
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_basic.RDS)."),
  make_option("--CVfiles", type = "character", help = "Path to the file containing predicted labels."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--method", type = "character", help = "Identifiers of the methods that were used")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


# Parameters
'%&%' = function(a,b) paste(a,b,sep="")
path = opt$outputDirec %&% opt$sampleName

################################################################################
## main code starts here
################################################################################
## load input data
IF_files <- list.files(opt$CVfiles,pattern = "Fold")
lfiles <- lapply(opt$method,function(x){IF_files[grep(x,IF_files)]})
r.files <- matrix(unlist(lfiles),nrow = length(lfiles[[1]]))
cv_files <- IF_files[grep(".RDS",IF_files)]
mat <- cbind(cv_files,r.files)
colnames(mat)[-1] <- opt$method

#load SCE
sce_data = readRDS(opt$SCE)
lab_data = colData(sce_data)$true_label
barcode = colnames(sce_data)

#data.frame: ground truth, predicted labels, CV fold
lab_mat <- as.data.frame(matrix(NA,
                                nrow = length(lab_data),
                                ncol = length(opt$method)+2))
lab_mat[,1] <- lab_data
for (i in 1:nrow(mat)){
  tmp <- readRDS(opt$CVfiles %&% mat[i,1])
  bc <- tmp$test_data
  idx <- which(barcode %in% bc)
  for (j in 1:length(opt$method)){
    tmp2 <- read_csv(opt$CVfiles %&% mat[i,j+1])
    lab_mat[idx,j+1] <- tmp2$x
  }
  lab_mat[idx,length(opt$method)+2] <- rep(i,length(idx))
}
colnames(lab_mat) <- c("lab_data",opt$method,"folds")
for (i in 1:length(opt$method)){
  lab_mat[,i+1] <- factor(lab_mat[,i+1],levels = levels(lab_mat[,1]))
}

#confusion matrix: accuracy, F1, overall statistics
accuracy <- matrix(NA,nrow = length(levels(lab_mat[,1])),
                   ncol = length(opt$method))
F1 <- matrix(NA,nrow = length(levels(lab_mat[,1])),
             ncol = length(opt$method))
overall <- matrix(NA,nrow = 7,
                  ncol = length(opt$method))
for (i in 1:length(opt$method)){
  assign(opt$method[i],confusionMatrix(lab_mat[,1],lab_mat[,i+1]))
  accuracy[,i] <- round(get(opt$method[i])$byClass[,"Balanced Accuracy"],4)
  F1[,i] <- round(get(opt$method[i])$byClass[,"F1"],4)
  overall[,i] <- round(get(opt$method[i])$overall,4)
}
colnames(accuracy) <- opt$method
rownames(accuracy) <- names(get(opt$method[1])$byClass[,"Balanced Accuracy"])
colnames(F1) <- opt$method
rownames(F1) <- names(get(opt$method[1])$byClass[,"F1"])
colnames(overall) <- opt$method
rownames(overall) <- names(get(opt$method[1])$overall)

#output
write.csv(accuracy,paste(path,"_byClass.csv",sep=""),quote = F,row.names = T)
write.csv(overall,paste(path,"_overall.csv",sep=""),quote = F,row.names = T)
write.csv(F1,paste(path,"_F1.csv",sep=""),quote = F,row.names = T)
