################################################################################
## Support Vector Machine Testing
################################################################################

#library
library(optparse)
library(tidyverse)
library(e1071)
library(cowplot)

#Data Path
#opt = list(
#  SCE = "/Users/bolars/Documents/celltyping/Zheng_sorted/ZhengSorted_data_test.RDS",
#  cell_labels = "/Users/bolars/Documents/celltyping/Zheng_sorted/ZhengSorted_lab_test.RDS",
#  outputDirec = "/Users/bolars/Documents/celltyping/Zheng_sorted/",
#  svm_model = "/Users/bolars/Documents/celltyping/Zheng_sorted/ZhengSorted_SVM_model.RDS",
#  sampleName = "ZhengSorted"
#)
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
  make_option("--svm_model", type = "character", help = "Path to the file containing the pretrained SVM model."),
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
sce_data = readRDS(opt$DataPath)
lab_data = readRDS(opt$LabelsPath)
svm_model = readRDS(opt$svm_model)
dat <- assay(sce_data)


#data.frame TO DO: updata to correct imput data format
xlevel = attr(svm_model$terms,"term.labels")
dat_filter <- dat[,colnames(dat) %in% xlevel]
data_svm <- bind_cols(lab_data,data_filter) %>%
  mutate(label=as.factor(label))

# SVM prediction
pred.svm <- predict(svm_model,newdata = data_svm,probability = T)
prob.svm <- attr(pred.svm,"probabilities")
idx.unk <- apply(prob.svm,1,max)<0.7
pred.svm <- as.character(pred.svm)
pred.svm[idx.unk] <- "unknown"
################################################################################
## save predicted labels
write.csv(pred.svm,path %&% "_SVM_pred_test_labels.csv")
