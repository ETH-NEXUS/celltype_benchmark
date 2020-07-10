################################################################################
## Fully Connected Neural Network Training
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
#  CVindex = "/Users/bolars/Documents/celltyping/benchmark_scripts/indexFold_01.RDS",
#  sampleName = "Fold_01",
#  method = "fcNN"
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
data_rf <- cbind(droplevels(lab_data),dat[,which(apply(dat,2,sum) != 0)])


#cross validation
training_fold = data_rf[rownames(data_rf) %in% cvindex$train_data, ]
test_fold = data_rf[rownames(data_rf) %in% cvindex$test_data, ]

x_train <- training_fold[,2:ncol(data_rf)]
y_train <- training_fold[,1]
x_test <- test_fold[,2:ncol(data_rf)]
y_test <- test_fold[,1]

to_one_hot <- function(labels, dimension = 6) {
  results <- matrix(0, nrow = length(labels), ncol = dimension)
  for (i in 1:length(labels))
    results[i, labels[[i]]] <- 1
  results
}
#shape
y_test <- as.matrix(as.numeric(y_test),ncol=1)
y_test <- to_one_hot(y_test)

y_train <- as.matrix(as.numeric(y_train),ncol=1)
y_train <- to_one_hot(y_train)

x_train <- as.matrix(x_train)
x_test <- as.matrix(x_test)

#fcNN
model <- keras_model_sequential() %>%
  layer_dense(units = 500, activation = 'relu', input_shape = c(ncol(data_rf))) %>%
  layer_dense(units = 200, activation = 'relu') %>%
  layer_dense(units = 50, activation = 'relu') %>%
  layer_dense(units = 6, activation = 'softmax')
#summary(model)
# compile model and intitialize weights
model %>% compile(
  optimizer = "adam",
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

if (!is_empty(cvindex$test_data)){
  history <- model %>% fit(
    x_train, y_train, 
    epochs = 20, batch_size=100,
    verbose=1, validation_data = list(x_test, y_test))
} else {
  history <- model %>% fit(
    x_train, y_train, 
    epochs = 20, batch_size=100,
    verbose=1)
}

#plot(history)
#results <- model %>% evaluate(x_test, y_test, verbose=0)
#results

if (!is_empty(cvindex$test_data)){
  tmp <- predict_classes(model,x_test)
  tmp2 <- tmp +1
  lev <- levels(test_fold[,1])
  res <- as.character(tmp2)
  for(i in 1:length(lev)){
    res[tmp2 %in% i] <- lev[i]
  }
  write.csv(res,path %&% "." %&% opt$method %&%
              "_predicted_label.csv",quote = F,row.names = F)
} else {
  save_model_hdf5(model, opt$outputDirec %&%
            opt$sampleName%&% "." %&%
            opt$method %&% "_model.h5")
}
