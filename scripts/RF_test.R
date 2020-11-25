################################################################################
## Random Forest Testing
################################################################################

#library
library(optparse)
library(tidyverse)
library(randomForest)
library(cowplot)
library(SingleCellExperiment)
#Data Path
#opt = list(
#  SCE = "/Users/bolars/Documents/celltyping/output/Newborn_PBMC_subset_test.RDS",
#  outputFile = "/Users/bolars/Documents/celltyping/output/training/pred_rf_newb_2.txt",
#  rf_model = "/Users/bolars/Documents/celltyping/output/training/All.RF_model_2.RDS",
#  sampleName = "RF_test",
#  method = "RF"
#)
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
  make_option("--rf_model", type = "character", help = "Path to the file containing the pretrained RF model."),
  make_option("--outputFile", type = "character", help = "Path to output file.")
  #make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name."),
  #make_option("--method", type = "character", help = "Method identifier. Attached to each output name.")
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
rf_model = readRDS(opt$rf_model)

#data frame
lab_data <- data.frame(label=lab_data)
dat <- as.data.frame(t(normcounts(sce_data)))
rownames(dat) <- colnames(sce_data)
colnames(dat) <- gsub("-","_",colnames(dat))
xlevel = names(rf_model$forest$ncat)
dat_filter <- dat[,colnames(dat) %in% xlevel]
missgene <- xlevel[!(xlevel %in% colnames(dat))]
if (!is_empty(missgene)){
  tmp <- matrix(0,nrow = nrow(dat_filter),ncol = length(missgene))
  colnames(tmp) <- missgene
  dat_filter <- cbind(dat_filter,tmp)
}
data_rf <- cbind(droplevels(lab_data),dat_filter)
data_rf$label <- as.factor(data_rf$label)

# Random Forest prediction
pred.rf <- predict(rf_model,newdata = data_rf[,-1])
prob.rf <- predict(rf_model,newdata = data_rf[,-1],type = "prob")
idx.unk <- apply(prob.rf,1,max)<0.3
pred.rf <- as.character(pred.rf)
pred.rf[idx.unk] <- "unknown"
pred.rf <- data.frame(celltype_final=pred.rf)

################################################################################
## save predicted labels
write_csv(pred.rf,opt$outputFile)
#write.csv(pred.rf,path %&% "." %&%
#            opt$method %&% "_predicted_test_labels.csv")
#confusionMatrix(data_rf[,1],factor(pred.rf,levels=levels(data_rf[,1])))
