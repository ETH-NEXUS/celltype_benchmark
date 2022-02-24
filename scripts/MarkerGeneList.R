################################################################################
## Finding marker genes for Garnett and scROSHI method
################################################################################

#library
library(optparse)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)

# command line arguments are parsed
option_list = list(
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--config", type = "character", help = "Path to cell type config file"),
  make_option("--format", type = "character", default = c("gmx","garnett"), help = "Output file format. (gmx, garentt)"),
  make_option("--foldChange", type = "double", default = c(2,1.5), help = "Limit testing to genes which showat least X-fold difference. (major type,sub type)"),
  make_option("--p_val_adj", type = "double", default = 0.05, help = "Threshold for adjusted p value")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

'%&%' = function(a,b) paste(a,b,sep="")

################################################################################
## main code starts here
################################################################################
## load input data
sce_data <- readRDS(opt$SCE)
#rownames(sce_data) <- rowData(sce_data)$SYMBOL #hotfix for HGNC symbols
lab_maj <- sce_data@metadata$ground_truth_major
lab_sub <- sce_data@metadata$ground_truth_minor

lab_sub[lab_sub=="none"] <- lab_maj[lab_sub=="none"]
lab_maj <- factor(lab_maj)
lab_sub <- factor(lab_sub)

#read config file
#################
conf <- read_tsv(opt$config, na = "none")
conf <- as.data.frame(conf)
subt <- vector(mode = "list",length = sum(!is.na(conf$Subtypes)))
names(subt) <- conf[which(!is.na(conf$Subtypes)),"Major_type"]
for (i in which(!is.na(conf$Subtypes))){
  subt[[which(names(subt)==conf[i,1])]] <- as.vector(str_split(conf[i,"Subtypes"],",",simplify = T))
}
majt <- conf$Major_type
#majt: major cell types (vector)
#subt: sub cell types (named list)

#checks if labels from sce data are also in the config file
stopifnot(all(levels(lab_maj) %in% majt),
          all(levels(lab_sub) %in% c(unlist(subt),majt))) 

#seurat data
seurat_data <- as.Seurat(sce_data, data = NULL)

#Major cell type labels; first run
Idents(seurat_data) = lab_maj
mgl <- vector(mode = "list",length = length(levels(lab_maj)))
names(mgl) <- levels(lab_maj)
for (grp in levels(lab_maj)){
  mark <- FindMarkers(seurat_data,
                      ident.1 = grp,
                      ident.2 = NULL,
                      only.pos = TRUE,
                      test.use = "wilcox",
                      logfc.threshold = log(opt$foldChange[1])) #hyperparameter 
  mgl[[grp]] <- rownames(mark[which(mark$p_val_adj<opt$p_val_adj),]) #hyperparameter
} # output mgl: list containing marker genes

#Sub cell types; second run
Idents(seurat_data) = lab_sub
#checks which of the major cell types (with subtypes) are present in the data
subtype <- subt[which(names(subt) %in% levels(lab_maj))]
#remove subtypes form list that are not present in the data
for(i in 1:length(subtype)){
  subtype[[i]] <- subtype[[i]][subtype[[i]] %in% levels(lab_sub)]
}
for (i in 1:length(subtype)){
  #subset with all subtypes of one major type
  seur_sub <- subset(seurat_data,idents = subtype[[i]])
  for (grp in subtype[[i]]){
    mark <- FindMarkers(seur_sub,
                        ident.1 = grp,
                        ident.2 = NULL,
                        only.pos = TRUE,
                        test.use = "wilcox",
                        logfc.threshold = log(opt$foldChange[2])) #hyperparameter
    mgl[[grp]] <- rownames(mark[which(mark$p_val_adj<opt$p_val_adj),]) #hyperparameter
  }
}
#gmx format
###########
if ("gmx" %in% opt$format){
  nRow <- max(as.numeric(summary(mgl)[,1]))+1
  mgl2 <- data.frame(matrix(NA,nrow = nRow,ncol = length(mgl)))
  for (i in 1:length(mgl)){
    mgl2[2:(length(mgl[[i]])+1),i] <- mgl[[i]]
  }
  colnames(mgl2) <- names(mgl)
  #output
  write_tsv(mgl2,opt$outputDirec %&% opt$sampleName %&% ".gmx",na="")
}

#garnett format
###############
if ("garnett" %in% opt$format){
  file.create(opt$outputDirec %&% opt$sampleName %&% ".garnett.tsv")
  #major types
  for (i in levels(lab_maj)){
    write_lines(paste(">",names(mgl[i]),sep=""),
                opt$outputDirec %&% opt$sampleName %&% ".garnett.tsv",
                append=T,sep="\n")
    write_lines(paste("expressed: ",paste(mgl[[i]],collapse=", "),sep=""),
                opt$outputDirec %&% opt$sampleName %&% ".garnett.tsv",
                append=T,sep="\n\n")
  }
  #subtypes
  for (i in unlist(subtype)){
    res <- lapply(subtype, function(ch) grep(names(mgl[i]), ch))
    nn <- names(subtype[sapply(res, function(x) length(x) > 0)])
    write_lines(paste(">",names(mgl[i]),sep=""),
                opt$outputDirec %&% opt$sampleName %&% ".garnett.tsv",
                append=T,sep="\n")
    write_lines(paste("expressed: ",paste(mgl[[i]],collapse=", "),sep=""),
                opt$outputDirec %&% opt$sampleName %&% ".garnett.tsv",
                append=T,sep="\n")
    write_lines(paste("subtype of: ",nn,sep=""),
                opt$outputDirec %&% opt$sampleName %&% ".garnett.tsv",
                append=T,sep="\n\n")
  }
}
