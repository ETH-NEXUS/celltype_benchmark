################################################################################
## Support Vector Machine Testing
################################################################################

#library
library(optparse)
library(tidyverse)
library(e1071)
library(cowplot)
library(SingleCellExperiment)
#Data Path
#opt = list(
#  SCE = "/Users/bolars/Documents/celltyping/benchmark_scripts/Zheng_merged_annotated.RDS",
#  outputFile = "/Users/bolars/Documents/celltyping/benchmark_scripts/pred_svm.txt",
#  svm_model = "/Users/bolars/Documents/celltyping/benchmark_scripts/All.SVM_model.RDS",
#  sampleName = "SVM_test",
#  method = "SVM"
#)
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
  make_option("--svm_model", type = "character", help = "Path to the file containing the pretrained SVM model."),
  make_option("--outputFile", type = "character", help = "Path to the output file.")
#  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name."),
#  make_option("--method", type = "character", help = "Method identifier. Attached to each output name.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Parameters
'%&%' = function(a,b) paste(a,b,sep="")
#path = opt$outputDirec %&% opt$sampleName

################################################################################
## main code starts here
################################################################################
## load input data
sce_data = readRDS(opt$SCE)
lab_maj <- sce_data@metadata$ground_truth_major
lab_data <- sce_data@metadata$ground_truth_minor
lab_data[lab_data=="none"] <- lab_maj[lab_data=="none"]
svm_model = readRDS(opt$svm_model)

#data frame
lab_data <- data.frame(label=lab_data)
dat <- as.data.frame(t(normcounts(sce_data)))
rownames(dat) <- colnames(sce_data)
colnames(dat) <- gsub("-","_",colnames(dat))
xlevel = attr(svm_model$terms,"term.labels")
dat_filter <- dat[,colnames(dat) %in% xlevel]
missgene <- xlevel[!(xlevel %in% colnames(dat))]
if (!is_empty(missgene)){
  tmp <- matrix(0,nrow = nrow(dat_filter),ncol = length(missgene))
  colnames(tmp) <- missgene
  dat_filter <- cbind(dat_filter,tmp)
}
data_svm <- cbind(droplevels(lab_data),dat_filter)
data_svm$label <- as.factor(data_svm$label)

# SVM prediction
pred.svm <- predict(svm_model,newdata = data_svm[,-1],probability = T)
prob.svm <- attr(pred.svm,"probabilities")
idx.unk <- apply(prob.svm,1,max)<0.5
pred.svm <- as.character(pred.svm)
pred.svm[idx.unk] <- "unknown"
################################################################################
## save predicted labels
write.csv(pred.svm,opt$outputFile)
#write.csv(pred.svm,path %&% "." %&%
#            opt$method %&% "_predicted_test_labels.csv")
#confusionMatrix(data_svm[,1],factor(pred.rf,levels=levels(data_svm[,1])))
