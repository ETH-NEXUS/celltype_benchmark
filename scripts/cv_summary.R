################################################################################
## Summary of the k-fold crossvalidation
################################################################################
#library
library(optparse)
library(tidyverse)
library(e1071)
library(cowplot)
library(caret)
library(Hmisc)
library(SingleCellExperiment)

#Data Path
# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_basic.RDS)."),
  make_option("--CVfiles", type = "character", help = "Path to the files containing predicted labels."),
  make_option("--CVindex", type = "character", help = "Path to the files containing cv index."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--method", type = "character", help = "Identifiers of the methods that were used")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
opt$method <- strsplit(opt$method, ",")[[1]]

# Parameters
'%&%' = function(a,b) paste(a,b,sep="")
path = opt$outputDirec %&% opt$sampleName

#functions
f.lab <- function(dat_lab,lab){
  #major labels
  lab <- unique(lab)
  notMaj <- which(!(dat_lab %in% lab[,1]))
  loop1 <- unique(dat_lab[notMaj])
  loop1 <- loop1[loop1!="Unknown"]
  res_mat <- data.frame(maj=dat_lab,sub=dat_lab)
  for(i in loop1){
    res_mat[res_mat[,1] == i,1] <- lab[lab[,2]==i,1]
  }
  #sub labels
  notSub <- which(!(dat_lab %in% lab[,2]))
  loop2 <- unique(dat_lab[notSub])
  loop2 <- loop2[loop2!="Unknown"]
  for(j in loop2){
    res_mat[res_mat[,2] == j,2] <- "Unknown"
  }
  return(res_mat)
}

################################################################################
## main code starts here
################################################################################

## load input data
IF_files <- list.files(opt$CVfiles,pattern = "Fold")
IF_files <- IF_files[grep(opt$sampleName,IF_files)]
IF_files <- IF_files[grep("predicted_label",IF_files)]
print(IF_files);print(rep("#",30))
lfiles <- lapply(opt$method,function(x){IF_files[grep(x,IF_files)]})
r.files <- matrix(unlist(lfiles),nrow = length(lfiles[[1]]))
cv_files <- list.files(opt$CVIndex,pattern = "Fold")
cv_files <- cv_files[grep(opt$sampleName,cv_files)]
print(cv_files)
mat <- cbind(cv_files,r.files)
colnames(mat)[-1] <- opt$method

#load SCE
sce_data = readRDS(opt$SCE)
lab_maj <- sce_data@metadata$ground_truth_major
lab_sub <- sce_data@metadata$ground_truth_minor
lab_sub[lab_sub=="none"] <- lab_maj[lab_sub=="none"]
barcode = colnames(sce_data)

#data.frame: ground truth, predicted labels, CV fold
lab_mat <- as.data.frame(matrix(NA,
                                nrow = length(lab_sub),
                                ncol = length(opt$method)+3))
lab_mat[,1] <- lab_maj
lab_mat[,2] <- lab_sub
for (i in 1:nrow(mat)){
  tmp <- readRDS(opt$CVIndex %&% mat[i,1])
  bc <- tmp$test_data
  idx <- which(barcode %in% bc)
  for (j in 1:length(opt$method)){
    tmp2 <- read_csv(opt$CVfiles %&% mat[i,j+1])
    lab_mat[idx,j+2] <- tmp2$x
  }
  lab_mat[idx,length(opt$method)+3] <- rep(i,length(idx))
}
colnames(lab_mat) <- c("lab_maj","lab_sub",opt$method,"folds")
print(head(lab_mat))

######################################
sub_lev <- unique(lab_mat$lab_sub)
maj_lev <- unique(lab_mat$lab_maj)
com_lev <- c(unique(c(lab_mat$lab_sub,lab_mat$lab_maj)),"Unknown")
accuracy_M <- matrix(NA,nrow = length(maj_lev),ncol = 0)
accuracy_S <- matrix(NA,nrow = length(sub_lev),ncol = 0)
F1_M <- matrix(NA,nrow = length(maj_lev),ncol = 0)
F1_S <- matrix(NA,nrow = length(sub_lev),ncol = 0)
overall_M <- matrix(NA,nrow = 8,ncol = 0)
overall_S <- matrix(NA,nrow = 8,ncol = 0)
nn <- character()

#confusion matrix: accuracy, F1, overall statistics
for (i in 1:length(opt$method)){
  j <- opt$method[i]
  if (all(lab_mat[,j] %in% sub_lev)){
    unkn <- 0
    pred <- factor(lab_mat[,j],levels = sub_lev)
    gt <- factor(lab_mat[,2],levels = sub_lev)
    t.tmp <- confusionMatrix(pred,gt)
    accuracy_S <- cbind(accuracy_S,round(t.tmp$byClass[,"Balanced Accuracy"],4))
    F1_S <- cbind(F1_S,round(t.tmp$byClass[,"F1"],4))
    over <- c(round(t.tmp$overall,4),unkn)
    names(over) <- c(names(t.tmp$overall),"Unknown")
    overall_S <- cbind(overall_S,over)
  }
  else if (all(lab_mat[,j] %in% com_lev)){
    res <- f.lab(lab_mat[,j],lab_mat[,1:2])
    #major
    t.dat <- res %>%
      add_column({{j}} :=lab_mat[,1])%>%
      filter(maj!="Unknown")
    unkn <- (nrow(res)-nrow(t.dat))/nrow(res)
    pred <- factor(t.dat[,1],levels = maj_lev)
    gt <- factor(t.dat[,j],levels = maj_lev)
    t.tmp <- confusionMatrix(pred,gt)
    accuracy_M <- cbind(accuracy_M,round(t.tmp$byClass[,"Balanced Accuracy"],4))
    F1_M <- cbind(F1_M,round(t.tmp$byClass[,"F1"],4))
    over <- c(round(t.tmp$overall,4),unkn)
    names(over) <- c(names(t.tmp$overall),"Unknown")
    overall_M <- cbind(overall_M,over)
    #sub
    t.dat <- res %>%
      add_column({{j}} :=lab_mat[,2])%>%
      filter(sub!="Unknown")
    unkn <- (nrow(res)-nrow(t.dat))/nrow(res)
    pred <- factor(t.dat[,2],levels = sub_lev)
    gt <- factor(t.dat[,j],levels = sub_lev)
    t.tmp <- confusionMatrix(pred,gt)
    accuracy_S <- cbind(accuracy_S,round(t.tmp$byClass[,"Balanced Accuracy"],4))
    F1_S <- cbind(F1_S,round(t.tmp$byClass[,"F1"],4))
    over <- c(round(t.tmp$overall,4),unkn)
    names(over) <- c(names(t.tmp$overall),"Unknown")
    overall_S <- cbind(overall_S,over)
    nn <- c(nn,j)
    }
}
#######################################
colnames(accuracy_S) <- opt$method
colnames(accuracy_M) <- nn
colnames(F1_S) <- opt$method
colnames(F1_M) <- nn
colnames(overall_S) <- opt$method
colnames(overall_M) <- nn

#output maj
write.csv(accuracy_M,paste(path,"_byClass_maj.csv",sep=""),quote = F,row.names = T)
write.csv(overall_M,paste(path,"_overall_maj.csv",sep=""),quote = F,row.names = T)
write.csv(F1_M,paste(path,"_F1_maj.csv",sep=""),quote = F,row.names = T)

#output maj
write.csv(accuracy_S,paste(path,"_byClass_sub.csv",sep=""),quote = F,row.names = T)
write.csv(overall_S,paste(path,"_overall_sub.csv",sep=""),quote = F,row.names = T)
write.csv(F1_S,paste(path,"_F1_sub.csv",sep=""),quote = F,row.names = T)
