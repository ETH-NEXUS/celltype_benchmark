################################################################################
## Finding marker genes for Garnet and NEXUS method
################################################################################

#library
library(optparse)
library(tidyverse)
library(Seurat)

#Data Path
opt = list(
  outputDirec = "/Users/bolars/Documents/celltyping/processed_sce/",
  SCE = "/Users/bolars/Documents/celltyping/processed_sce/Zheng_sorted_merged.genes_cells_filtered.corrected.ground-truth.RDS",
  sampleName = "MakerGeneList"
)
# command line arguments are parsed
#option_list = list(
#  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
#  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
#  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name.")
#)

#opt_parser = OptionParser(option_list = option_list)
#opt = parse_args(opt_parser)
'%&%' = function(a,b) paste(a,b,sep="")

################################################################################
## main code starts here
################################################################################
## load input data
sce_data <- readRDS(opt$SCE)
lab_data <- colData(sce_data)$true_label

seutat_data <- as.Seurat(sce_data, data = NULL)
Idents(seutat_data) = lab_data

mgl <- vector(mode = "list",length = length(levels(lab_data)))
names(mgl) <- levels(lab_data)

for (grp in levels(lab_data)){
  mark <- FindMarkers(seutat_data,
                      ident.1 = grp,
                      ident.2 = NULL,
                      only.pos = TRUE,
                      test.use = "wilcox",
                      logfc.threshold = 1,
                      min.cells.feature = 20,
                      min.cells.group = 20)
  mgl[[grp]] <- rownames(mark[which(mark$p_val_adj<0.05),])
}

#gmx format
nRow <- max(as.numeric(summary(mgl)[,1]))+1
mgl2 <- data.frame(matrix(NA,nrow = nRow,ncol = length(mgl)))
for (i in 1:length(mgl)){
  mgl2[2:(length(mgl[[i]])+1),i] <- mgl[[i]]
}
colnames(mgl2) <- names(mgl)

#output
write_tsv(mgl2,opt$outputDirec %&% opt$sampleName %&% ".gmx",na="")


