################################################################################
## Random Forest Training
################################################################################

#library
library(optparse)
library(tidyverse)
library(randomForest)
library(cowplot)
library(caret)


#Data Path
opt = list(
  SCE = "/Users/bolars/Documents/celltyping/Zheng_sorted/ZhengSorted_data.RDS",
  cell_labels = "/Users/bolars/Documents/celltyping/Zheng_sorted/ZhengSorted_lab.RDS",
  outputDirec = "/Users/bolars/Documents/celltyping/benchmark_scripts/",
  CVindex = "/Users/bolars/Documents/celltyping/benchmark_scripts/indexAll.RDS",
  sampleName = "indexAll_RF"
)
# command line arguments are parsed
#option_list = list(
#  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
#  make_option("--cell_labels", type = "character", help = "Path to the file containing the cell labels."),
#  make_option("--CVindex", type = "character", help = "Path to the file containing the indices of the cross validation."),
#  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
#  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name.")
#)

#opt_parser = OptionParser(option_list = option_list)
#opt = parse_args(opt_parser)


# Parameters
'%&%' = function(a,b) paste(a,b,sep="")
path = opt$outputDirec %&% opt$sampleName
Ntrees <- 500

################################################################################
## main code starts here
################################################################################
## load input data
sce_data = readRDS(opt$SCE)
lab_data = readRDS(opt$cell_labels)
cvindex = readRDS(opt$CVindex)
#dat <- assay(sce_data)

#data frame TO DO: updata to correct imput data format
data_rf <- bind_cols(lab_data,sce_data) %>%
  mutate(label=as.factor(label))
lev = levels(data_rf$label)

#cross validation
training_fold = data_rf[cvindex$train_data, ]
test_fold = data_rf[cvindex$test_data, ]

#random forest
classifier <- randomForest(label~.,
                           data=training_fold,
                           ntree=Ntrees)
if (!is_empty(cvindex$test_data)){
  y_pred = predict(classifier, newdata = test_fold[,-1])
  #y_pred <- factor(y_pred,levels = lev)
  #y_test <- test_fold[, 1] %>%
  #mutate(label=factor(label,levels = lev))
  #confusion matrix
  #cm = confusionMatrix(y_test$label, y_pred)
  #saveRDS(cm, path %&% "_cm.RDS")
  write.csv(y_pred,path %&% "_pred_label.csv",quote = F,row.names = F)
} else {
  saveRDS(classifier, path %&% "_model.RDS")
}