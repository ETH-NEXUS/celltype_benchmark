#add_sce_groundtruth <- function(rds_path, labels_path){
library("future")

sce <- readRDS(rds_path)
sce@metadata$ground_truth_major <- runif(ncol(sce))
sce@metadata$ground_truth_minor <- runif(ncol(sce))

labels <- read.table(labels_path, skip=1, sep ='\t')
labels <- labels[labels$V2!="unannotated",]
labels <- labels[labels$V2!="B cell T cell doublet",]

labels_full <- labels$V2
labels_short <- labels$V3
names(labels_full) <- labels$V1
names(labels_short) <- labels$V1

run_time_test <- function(tool){
barcode_index <- c()
time_start <- Sys.time()
tool()
time_end <- Sys.time()
time_elapsed <- time_start - time_end
print(time_elapsed)
}
    
match_sequential <- function(){	#time elapsed: 4 mins (now 3 mins)
plan(future::sequential)
i <- 1
i_max <- length(names(labels_full))
num_matched <- 0
num_unmatched <- 0
barcode_index <- vector("list", i_max)
    
    for(barcode in names(labels_full)){
    barcode_trunc <- substr(barcode,1,16)
    if(barcode_trunc %in% sce$barcodes){
        #print(paste(barcode, as.character(labels_full[barcode]), sep = " --- "))
        match_index <- which(barcode_trunc == sce$barcodes)
        names(match_index) <- barcode_trunc
        barcode_index[[i]] <- match_index
        #print(paste(which(barcode == names(labels_full)), length(names(labels_full)), sep = "/"))
        num_matched <- num_matched + 1
    } else {
        num_unmatched <- num_unmatched + 1}
    i <- i + 1
}
    print(num_matched)
	print(num_unmatched)
    return(unlist(barcode_index))
	}

match_parallel_lapply <- function(){	#time elapsed: 1.5 mins (now 5 mins)
plan(future::multiprocess)
i <- 1
i_max <- length(names(labels_full))
num_matched <- 0
num_unmatched <- 0
barcode_index <- vector("list", i_max)
	
future.apply::future_lapply(names(labels_full), function(barcode){
    barcode_trunc <- substr(barcode,1,16)
    if(barcode_trunc %in% sce$barcodes){
        #print(paste(which(barcode == names(labels_full)), length(names(labels_full)), sep = "/"))
        match_index <- which(barcode_trunc == sce$barcodes)
        names(match_index) <- barcode
        barcode_index[[i]] <- match_index
        num_matched <- num_matched + 1
    } else {
        num_unmatched <- num_unmatched + 1}
    i <- i + 1
    })
    print(num_matched)
	print(num_unmatched)
    }

# Other approach:
#smallFrame$newColumn <- unlist(someList[match(smallFrame$colA, names(someList))])
match_index <- match(sce$barcodes, substr(names(labels_full),1,16)) #SIGNIFICANTLY FASTER
#equivalent to: x %in% table except returns index

sce@metadata$ground_truth_major <- as.character(labels_short[match_index])
sce@metadata$ground_truth_minor <- as.character(labels_full[match_index])

#bad matches/NA; ignored commented chunk below this
which(is.na(as.character(labels_full[match_index])))
#no issue when done like this: e.g. as.character(labels_full[116])


#REAL ISSUE index has numbers larger than actual list - why?
# REAL REAL ISSUE: I CHANGED LABELS LIST TO EXCLUDE CELL TYPES WITHOUT REDOING MATCHING!!!)

#na_index <- is.na(sce@metadata$ground_truth_major)
#print(paste(sum(na_index), "cells removed from SCE."))
#na_barcodes <- sce@barcodes[na_index]
#row.names(labels) <- substr(labels$V1,1,16)
#print(table(labels[na_barcodes,]$V2))
#}

#Next issue: in SCE, cell types seem to be arranged sequentially; would not make sense since barcodes are not sorted by cell type...or are they...?
#However, they seem to match the list barcodes & labels. Maybe the subsetting of data by cell type caused barcodes to structure in this way & post-processing doesn't arrange them.


#sce[,barcode_trunc]@metadata$ground_truth_major <- as.character(labels_full[barcode])
#sce[,barcode_trunc]@metadata$ground_truth_minor <- as.character(labels_short[barcode])
    
#sce@metadata$ground_truth_major != NULL 
    }

dir_base <- "/cluster/project/nexus/marcus/cell_typing_benchmark/pbmc_label_lists/"
adult_labels <- paste(dir_base, "2020-Mar-Census-Adult-Immune-10x.ground_truth.tsv", sep = "")
adult_sce <- "/cluster/project/nexus/benchmarking/celltyping/data/subset/2020-Mar-Census-Adult-Immune-10x.scp.genes_cells_filtered.corrected.RDS"
#zheng_labels <- paste(dir_base, "Zheng_sorted_merged.ground_truth.tsv", sep = "")

#newborn_labels <- paste(dir_base, "2020-Mar-Census-Newborn-Blood-10x.ground_truth.tsv", sep = "")

add_sce_groundtruth(adult_sce, adult_labels)
