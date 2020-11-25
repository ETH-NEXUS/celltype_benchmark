################################################################################
## Support Vector Machine Training
################################################################################

#library
library(optparse)
library(tidyverse)
library(e1071)
library(cowplot)
library(caret)
library(SingleCellExperiment)

#Data Path
#opt = list(
#  SCE = "/Users/bolars/Documents/celltyping/output/Adult_PBMC_subset_test.RDS",
#  outputDirec = "/Users/bolars/Documents/celltyping/output/training/",
#  CVindex = "/Users/bolars/Documents/celltyping//output/training/Adult_PBMC_subset.index_All.RDS",
#  sampleName = "All",
#  method = "SVM"
#)
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_basic.RDS)."),
  make_option("--CVindex", type = "character", help = "Path to the file containing the indices of the cross validation."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--method", type = "character", help = "Method identifier. Attached to each output name.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Parameters
'%&%' = function(a,b) paste(a,b,sep="")

################################################################################
## main code starts here
################################################################################
## load input data
sce_data = readRDS(opt$SCE)
lab_maj <- sce_data@metadata$ground_truth_major
lab_data <- sce_data@metadata$ground_truth_minor
lab_data[lab_data=="none"] <- lab_maj[lab_data=="none"]
cvindex = readRDS(opt$CVindex)

#data frame
lab_data <- data.frame(label=lab_data)
dat <- as.data.frame(t(normcounts(sce_data)))
rownames(dat) <- colnames(sce_data)
colnames(dat) <- gsub("-","_",colnames(dat))
data_svm <- cbind(droplevels(lab_data),dat[,which(apply(dat,2,sum) != 0)])
data_svm$label <- as.factor(data_svm$label)

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
  write.csv(y_pred,opt$outputDirec %&%
              opt$sampleName %&% "." %&%
              opt$method %&% "_predicted_labels.csv",
            quote = F,row.names = F)
} else {
  saveRDS(classifier, opt$outputDirec %&%
            opt$sampleName %&% "." %&%
            opt$method %&% "_model.RDS")
}
