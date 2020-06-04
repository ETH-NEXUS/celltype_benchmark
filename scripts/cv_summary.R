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
#  sampleName = "CV_summary"
#)
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_basic.RDS)."),
  make_option("--CVfiles", type = "character", help = "Path to the file containing predicted labels."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name.")
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
IF_files <- list.files(opt$CVfiles,pattern = "indexFold")
rf_files <- IF_files[grep("_RF_",IF_files)]
svm_files <- IF_files[grep("_SVM_",IF_files)]
cv_files <- IF_files[grep(".RDS",IF_files)]
mat <- cbind(cv_files,rf_files,svm_files)

#load data
sce_data = readRDS(opt$SCE)
lab_data = colData(sce_data)$true_label
barcode = colnames(sce_data)

lab_mat <- data.frame(true_lab=lab_data,rf_lab=NA,svm_lab=NA,folds=NA)
for (i in 1:nrow(mat)){
  tmp <- readRDS(opt$CVfiles %&% mat[i,1])
  bc <- tmp$test_data
  idx <- which(barcode %in% bc)
  rf <- read_csv(opt$CVfiles %&% mat[i,2])
  svm <- read_csv(opt$CVfiles %&% mat[i,3])
  lab_mat[idx,"rf_lab"] <- rf$x
  lab_mat[idx,"svm_lab"] <- svm$x
  lab_mat[idx,"folds"] <- rep(i,length(idx))
}
lab_mat<- lab_mat %>%
  mutate(true_lab=factor(true_lab),
         rf_lab=factor(rf_lab),
         svm_lab=factor(svm_lab))
cm_rf <- confusionMatrix(lab_mat$true_lab,lab_mat$rf_lab)
cm_svm <- confusionMatrix(lab_mat$true_lab,lab_mat$svm_lab)

accuracy <- cbind(round(cm_rf$byClass[,"Balanced Accuracy"],4),
                  round(cm_svm$byClass[,"Balanced Accuracy"],4))
F1 <- cbind(round(cm_rf$byClass[,"F1"],4),
                  round(cm_svm$byClass[,"F1"],4))
overall <- cbind(round(cm_rf$overall,4),
                 round(cm_svm$overall,4))
colnames(accuracy) <- c("RF","SVM")
colnames(overall) <- c("RF","SVM")
colnames(F1) <- c("RF","SVM")
write.csv(accuracy,paste(path,"_byClass.csv",sep=""),quote = F,row.names = T)
write.csv(overall,paste(path,"_overall.csv",sep=""),quote = F,row.names = T)
write.csv(F1,paste(path,"_F1.csv",sep=""),quote = F,row.names = T)