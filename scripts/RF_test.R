################################################################################
## Random Forest Testing
################################################################################

#library
library(optparse)
library(tidyverse)
library(randomForest)
library(cowplot)

#Data Path
opt = list(
  DataPath = "/Users/bolars/Documents/celltyping/Zheng_sorted/ZhengSorted_data_test.RDS",
  LabelsPath = "/Users/bolars/Documents/celltyping/Zheng_sorted/ZhengSorted_lab_test.RDS",
  outputDirec = "/Users/bolars/Documents/celltyping/Zheng_sorted/",
  rf_model = "/Users/bolars/Documents/celltyping/Zheng_sorted/ZhengSorted_RF_model.RDS",
  sampleName = "ZhengSorted_RF"
)
# command line arguments are parsed
#option_list = list(
#  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
#  make_option("--cell_labels", type = "character", help = "Path to the file containing the cell labels."),
#  make_option("--rf_model", type = "character", help = "Path to the file containing the pretrained RF model."),
#  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
#  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name.")
#)

#opt_parser = OptionParser(option_list = option_list)
#opt = parse_args(opt_parser)


# Parameters
'%&%' = function(a,b) paste(a,b,sep="")
path = opt$outputDirec %&% opt$sampleName

################################################################################
## main code starts here
################################################################################
## load input data
sce_data = readRDS(opt$DataPath)
lab_data = readRDS(opt$LabelsPath)
rf_model = readRDS(opt$rf_model)
dat <- assay(sce_data)

#data.frame TO DO: updata to correct imput data format
xlevel = names(rf_model$forest$ncat)
dat_filter <- dat[,colnames(dat) %in% xlevel]
data_rf <- bind_cols(lab_data,data_filter) %>%
  mutate(label=as.factor(label))

# Random Forest prediction
pred.rf <- predict(rf_model,newdata = data_rf)
prob.rf <- predict(rf_model,newdata = data_rf,type = "prob")
idx.unk <- apply(prob.rf,1,max)<0.7
pred.rf <- as.character(pred.rf)
pred.rf[idx.unk] <- "unknown"

################################################################################
## save predicted labels
write.csv(pred.rf,path %&% "_pred_test_labels.csv")