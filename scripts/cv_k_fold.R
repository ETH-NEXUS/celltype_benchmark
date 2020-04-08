################################################################################
## creating the k-fold crossvalidation data
################################################################################

#library
library(optparse)
library(tidyverse)
library(caret)

# command line arguments are parsed
option_list = list(
  make_option("--cell_labels", type = "character", help = "Path to the file containing cell lables."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Parameters
Kfold <- 10

################################################################################
## main code starts here
################################################################################
## load input data
lab_data <- read_csv(opt$LabelsPath)
lab_data <- rename(lab_data,"label"="x")

# cross validation
##################
folds = createFolds(lab_data$label, k = Kfold,list = T)
allIndex <- 1:nrow(lab_data)
for(i in 1:length(folds)){
  test_data <- folds[[i]]
  train_data <- allIndex[-folds[[i]]]
  saveRDS(list(test_data=test_data,train_data=train_data),
            paste(opt$outputDirec,"index",names(folds)[i],".RDS",sep=""))
}
saveRDS(list(test_data=numeric(),train_data=allIndex),paste(opt$outputDirec,"indexAll.RDS",sep=""))
