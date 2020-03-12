##### Test run for CellAssign #####
reticulate::use_condaenv("cellassign2",required=TRUE)
library(cellassign)
library(Matrix)
library(SingleCellExperiment)
library(org.Hs.eg.db)

#matrix_dir = "/cluster/project/nexus/marcus/cell_typing_benchmark/Cross_Validation/filtered_feature_bc_matrix/"
matrix_dir = "/cluster/project/nexus/benchmarking/celltyping/test_pilot/cellranger_run/b_cells/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv") #or .tsv.gz; make a robust fxn
features.path <- paste0(matrix_dir, "features.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

sce <- SingleCellExperiment(mat)
names(sce@assays) <- "counts"

marker_path <- "/cluster/project/nexus/benchmarking/celltyping/marker_files/PBMC/PBMC_markers.txt" #don't use CellAssign specific marker file...needs reformatting
marker_file <-t(read.csv(marker_path, header=FALSE, row.names=1))

# Linda's PBMC marker file - made manually first
marker_list <- list(
  "B.cells.naive_normal_Newman15" = c(as.character(marker_file[,1])),                                    
  "B.cells.memory_normal_Newman15" = c(as.character(marker_file[,2])),                                   
  "Plasma.cells_normal_Newman15" = c(as.character(marker_file[,3])),                                      
  "T.cells.CD8_normal_Newman15" = c(as.character(marker_file[,4])),                                       
  "T.cells.CD4.naive_normal_Newman15" = c(as.character(marker_file[,5])),                                 
  "T.cells.CD4.memory.resting_normal_Newman15" = c(as.character(marker_file[,6])),                        
  "T.cells.CD4.memory.activated_normal_Newman15" = c(as.character(marker_file[,7])),                      
  "T.cells.follicular.helper_normal_Newman15" = c(as.character(marker_file[,8])),                         
  "T.cells.regulatory_normal_Newman15" = c(as.character(marker_file[,9])),                                
  "T.cells.gamma.delta_normal_Newman15" = c(as.character(marker_file[,10])),                               
  "NK.cells.resting_normal_Newman15" = c(as.character(marker_file[,11])),                                  
  "NK.cells.activated_normal_Newman15" = c(as.character(marker_file[,12])),                                
  "Monocytes_normal_Newman15" = c(as.character(marker_file[,13])),                                         
  "Macrophages.M0_normal_Newman15" = c(as.character(marker_file[,14])),                                    
  "Macrophages.M1_normal_Newman15" = c(as.character(marker_file[,15])),                                    
  "Macrophages.M2_normal_Newman15" = c(as.character(marker_file[,16])),                                    
  "Dendritic.cells.resting_normal_Newman15" = c(as.character(marker_file[,17])),                           
  "Dendritic.cells.activated_normal_Newman15" = c(as.character(marker_file[,18])),                         
  "Mast.cells.resting_normal_Newman15" = c(as.character(marker_file[,19])),                                
  "Mast.cells.activated_normal_Newman15" = c(as.character(marker_file[,20])),                              
  "Eosinophils_normal_Newman15" = c(as.character(marker_file[,21])),                                       
  "Neutrophils_normal_Newman15" = c(as.character(marker_file[,22])),                                       
  "Dendritic.cells_normal_Newman15.MERGED" = c(as.character(marker_file[,23])),                            
  "T.cells_melanoma_Tirosh16" = c(as.character(marker_file[,24])),                                         
  "B.cells_melanoma_Tirosh16" = c(as.character(marker_file[,25])),                                         
  "Macrophages_melanoma_Tirosh16" = c(as.character(marker_file[,26])),                                     
  "Endothelial.cells_melanoma_Tirosh16" = c(as.character(marker_file[,27])),                               
  "Naive.CD8pos.T.cell.Peripheral.blood_Normal_Zhang18" = c(as.character(marker_file[,28])),               
  "Effector.CD8+.memory.T.(Tem).cell.Peripheral.blood_Normal_Zhang18" = c(as.character(marker_file[,29])), 
  "Naive.CD4+.T.cell.Peripheral.blood_Normal_Zhang18" = c(as.character(marker_file[,30])),                 
  "Regulatory.T.(Treg).cell.Peripheral.blood_Normal_Zhang18" = c(as.character(marker_file[,31])),          
  "CD141+CLEC9A+.dendritic.cell.Blood_Normal_Zhang18" = c(as.character(marker_file[,32])),                 
  "CD1C+_A.dendritic.cell.Blood_Normal_Zhang18" = c(as.character(marker_file[,33])),                       
  "CD1C+_B.dendritic.cell.Blood_Normal_Zhang18" = c(as.character(marker_file[,34])),                       
  "CD1C-CD141-.dendritic.cell.Blood_Normal_Zhang18" = c(as.character(marker_file[,35])),                   
  "AXL+SIGLEC6+.dendritic.cell.Blood_Normal_Zhang18" = c(as.character(marker_file[,36])),                  
  "Plasmacytoid.dendritic.cell.Blood_Normal_Zhang18" = c(as.character(marker_file[,37])) 
)

marker_mat <- marker_list_to_mat(marker_list)
gene_symbols <- as.vector(mapIds(org.Hs.eg.db, keys = row.names(marker_mat), keytype = "SYMBOL", column="ENSEMBL"))
row.names(marker_mat) <- gene_symbols
marker_mat <- marker_mat[row.names(marker_mat) %in% row.names(sce),]

pbmc_cds <- sce[row.names(marker_mat),]

valid_cells <- c() # valid cells: cells w/ non-zero counts for marker genes
for(cell in colnames(pbmc_cds)){
  if(sum(counts(pbmc_cds)[,cell]) > 0){
    valid_cells <- c(valid_cells,cell)
  }
}

pbmc_cds <- pbmc_cds[, valid_cells]

s <- sizeFactors(sce)

# It is critical that gene expression data containing only marker genes is used as input to cellassign.
fit <- cellassign(exprs_obj = pbmc_cds, 
                  marker_gene_info = marker_mat, 
                  s = s, 
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)

print(fit)
print(head(celltypes(fit)))
print(str(mleparams(fit)))
