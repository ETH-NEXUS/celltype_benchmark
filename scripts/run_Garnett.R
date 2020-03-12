args <- commandArgs(TRUE)

run_Garnett_CV <- function(DataPath, LabelsPath, CV_RDataPath, GenesPath, MarkerPath, OutputDir, Human){
  "
  run Garnett
  Wrapper script to run Garnett on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  GenesPath : Path to the file with the genenames
  MarkerPath : Path to the file with marker genes
  OutputDir : Output directory defining the path of the exported file.
  Human : boolean indicating whether the dataset is human (TRUE) or mouse (FALSE)
  "


library(garnett)
library(org.Hs.eg.db)

#data_path <- "/cluster/project/nexus/benchmarking/celltyping/test_pilot/cellranger_run/b_cells"
#marker_path <- "/cluster/project/nexus/benchmarking/celltyping/garnett_files/marker_files/hsPBMC_markers.txt"

pbmc_cds <- load_cellranger_data(data_path)
barcodes <- as.vector(pData(pbmc_cds)$barcode)

process_cds_monocle <- function(cds){
  cds <- preprocess_cds(cds, num_dim = 100)
  cds <- reduce_dimension(cds)
  return(cds)
}

pbmc_cds <- process_cds_monocle(pbmc_cds)

training_indices <- sample(c(TRUE,FALSE), length(barcodes), replace=TRUE, prob=c(0.7,0.3))
train_set <- barcodes[training_indices]
test_set <- barcodes[!training_indices]

pbmc_train <- pbmc_cds[,train_set]
pbmc_train <- process_cds_monocle(pbmc_train)
pbmc_test <- pbmc_cds[,test_set]
pbmc_test <- process_cds_monocle(pbmc_test)

pbmc_classifier <- train_cell_classifier(cds = pbmc_train, 
                                         marker_file = marker_path,
                                         db=org.Hs.eg.db,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")


pbmc_test <- classify_cells(pbmc_test, 
                            pbmc_classifier, 
                            db = org.Hs.eg.db, 
                            cluster_extend = TRUE,
                            cds_gene_id_type = "ENSEMBL")


#true_labels <- list(lab_test)
pred_labels <- list(pData(pbmc_test)$cluster_ext_type)

for (i in c(1:n_folds)){
  lab_train = labels[Train_Idx[[i]]]
  lab_test = labels[Test_Idx[[i]]]
  
  train = data[,Train_Idx[[i]]]
  test = data[,Test_Idx[[i]]]
  
  cells_train = cells[Train_Idx[[i]]]
  cells_test = cells[Test_Idx[[i]]]
  
  pdata_train = data.frame(cells_train)
  pdata_test = data.frame(cells_test)
  
  row.names(train) <- row.names(fdata)
  row.names(test) <- row.names(fdata)
  colnames(train) <- row.names(pdata_train)
  colnames(test) <- row.names(pdata_test)
  
  #pd_train <- new("data.frame", data = pdata_train)
  #pd_test <- new("data.frame", data = pdata_test)
  pd_train <- pdata_train
  pd_test <- pdata_test
  "Generating CDS objects..."
  pbmc_cds_train <- new_cell_data_set(as(train, "dgCMatrix"), cell_metadata = pd_train, gene_metadata = fd)
  pbmc_cds_test <- new_cell_data_set(as(test, "dgCMatrix"), cell_metadata = pd_test, gene_metadata = fd)
  
  pbmc_cds_train <- estimate_size_factors(pbmc_cds_train)
  pbmc_cds_test <- estimate_size_factors(pbmc_cds_test)
  
  # training
  print("Starting training...")
  start_train <- Sys.time()
  
  pbmc_classifier <- train_cell_classifier(cds = pbmc_cds_train, 
                                           marker_file = MarkerPath,
                                           db=org.Hs.eg.db,
                                           cds_gene_id_type = "SYMBOL",
                                           num_unknown = 50,
                                           marker_file_gene_id_type = "SYMBOL")
  end_train <- Sys.time()
  train_time[i] <- as.numeric(end_train - start_train)
  
  # testing
  print("Starting testing...")
  start_test <- Sys.time()
  
  pbmc_cds_test <- classify_cells(pbmc_cds_test, 
                                  pbmc_classifier, 
                                  db = org.Hs.eg.db, 
                                  cluster_extend = TRUE,
                                  cds_gene_id_type = "SYMBOL")
  end_test <- Sys.time()
  test_time[i] <- as.numeric(end_test - start_test)
  
  true_labels[i] <- list(lab_test)
  pred_labels[i] <- list(pData(pbmc_cds_test)$cluster_ext_type)
  print(paste("Fold",i,"completed.",sep = " "))
  
}
}

run_Garnett_CV(args[1], args[2], args[3], args[4], args[5], args[6], args[7])

