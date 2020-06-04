################################################################################
## Support Vector Machine Training
################################################################################

#library
library(optparse)
library(tidyverse)
library(e1071)
library(cowplot)
library(caret)


#Data Path
#opt = list(
#  SCE = "/Users/bolars/Documents/celltyping/benchmark_scripts/Zheng_sorted_merged.genes_cells_filtered.corrected.ground-truth.RDS",
#  outputDirec = "/Users/bolars/Documents/celltyping/benchmark_scripts/",
#  CVindex = "/Users/bolars/Documents/celltyping/benchmark_scripts/indexFold01.RDS",
#  sampleName = "Fold01"
#)
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_basic.RDS)."),
  make_option("--CVindex", type = "character", help = "Path to the file containing the indices of the cross validation."),
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
sce_data = readRDS(opt$SCE)
lab_data = colData(sce_data)$true_label
cvindex = readRDS(opt$CVindex)

#data frame
lab_data <- data.frame(label=lab_data)
dat <- as.data.frame(t(normcounts(sce_data)))
rownames(dat) <- colnames(sce_data)
colnames(dat) <- gsub("-","_",colnames(dat))
data_svm <- cbind(droplevels(lab_data),dat[,which(apply(dat,2,sum) != 0)])

#cross validation
training_fold = data_svm[rownames(data_svm) %in% cvindex$train_data, ]
test_fold = data_svm[rownames(data_svm) %in% cvindex$test_data, ]

#SVM
classifier <- svm(label ~ .,
                 data = training_fold,
                 method = "C-classification",
                 kernel = "linear",
                 probability=T)
if (!is_empty(cvindex$test_data)){
  y_pred = predict(classifier, newdata = test_fold[,-1])
  #y_pred <- factor(y_pred,levels = lev)
  #y_test <- test_fold[, 1] %>%
    #mutate(label=factor(label,levels = lev))
  #confusion matrix
  #cm = confusionMatrix(y_test$label, y_pred)
  #saveRDS(cm, path %&% "_cm.RDS")
  write.csv(y_pred,path %&% "_SVM_pred_label.csv",quote = F,row.names = F)
} else {
  saveRDS(classifier, path %&% "_SVM_model.RDS")
}