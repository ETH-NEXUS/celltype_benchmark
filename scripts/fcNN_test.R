################################################################################
## Fully Connected Neural Network Testing
################################################################################

#library
library(optparse)
library(tidyverse)
library(caret)
library(cowplot)
library(keras)

#Data Path
#opt = list(
#  SCE = "/Users/bolars/Documents/celltyping/benchmark_scripts/Zheng_sorted_merged.genes_cells_filtered.corrected.ground-truth.RDS",
#  outputDirec = "/Users/bolars/Documents/celltyping/benchmark_scripts/",
#  nn_model = "/Users/bolars/Documents/celltyping/benchmark_scripts/All.fcNN_model.h5",
#  sampleName = "fcNN_Test",
#  method = "fcNN"
#)
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
  make_option("--nn_model", type = "character", help = "Path to the file containing the pretrained fcNN model."),
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
lab_data = colData(sce_data)$true_label
nn_model = load_model_hdf5(opt$nn_model)
dat <- as.data.frame(t(normcounts(sce_data)))
rownames(dat) <- colnames(sce_data)
colnames(dat) <- gsub("-","_",colnames(dat))
data_nn <- cbind(droplevels(lab_data),dat[,which(apply(dat,2,sum) != 0)])

x_test <- data_nn[,2:ncol(data_nn)]
x_test <- as.matrix(x_test)

#predict labels
tmp <- predict_classes(nn_model,x_test)
tmp2 <- tmp +1
lev <- levels(test_fold[,1])
res <- as.character(tmp2)
for(i in 1:length(lev)){
  res[tmp2 %in% i] <- lev[i]
}
write.csv(res,opt$outputFile)
#write.csv(res,path %&% "." %&% opt$method %&%
#            "_predicted_test_label.csv",quote = F,row.names = F)
#confusionMatrix(data_nn[,1],factor(res,levels=levels(data_nn[,1])))
