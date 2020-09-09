# Subset counts data to 10% of original size, respecting proportions of cell types featured in a ground-truth labels file.

library(optparse)

option_list = list(
  make_option("--input_data", type = "character", help = "Path to hdf5 file with raw counts data."),
  make_option("--true_labels", type = "character", help = "Path to the file containing the ground-truth labels with cell types."),
  make_option("--output_subset", type = "character", help = "Path to the directory where output files will be written.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

###

select_barcodes <- function(labels_path){
    
    labels <- read.table(labels_path, skip=1, sep ='\t')
    final_barcodes <- c()
    
    proportions <- table(labels$V2)
        for(cell_type in names(proportions)){
            cell_count <- proportions[[cell_type]]
            if(cell_count > 0){
                cell_select <- cell_count / 10
                barcodes <- labels[labels$V2==cell_type,]$V1
                subset_barcodes <- as.character(sample(barcodes, cell_select))
                subset_barcodes <- substr(subset_barcodes,1,16) #only for hdf5 files
                final_barcodes <- c(final_barcodes, subset_barcodes)
                print(head(final_barcodes))
                }       
    		}
    return(final_barcodes)
}

barcodes_subset <- select_barcodes(opt$true_labels)

##################################################

create_hdf5_subset <- function(barcodes_subset, hdf5_in, hdf5_out){
    
    cells <- h5read(hdf5_in, "cell_attrs/cell_names")
    cells_index <- match(barcodes_subset, cells) #none are matching...?
	h5createFile(hdf5_out)
	h5createGroup(hdf5_out, "cell_attrs")
	h5createGroup(hdf5_out, "gene_attrs")

	h5write(as.character(cells[cells_index]), file = hdf5_out, "cell_attrs/cell_names")
	h5write(as.character(h5read(hdf5_in, "cell_attrs/cells_on_rows")), file = hdf5_out, "cell_attrs/cells_on_rows") #can I remove this? issue w/ factors
	h5write(h5read(hdf5_in, "gene_attrs/gene_ids"), file = hdf5_out, "gene_attrs/gene_ids")
	h5write(h5read(hdf5_in, "gene_attrs/gene_names"), file = hdf5_out, "gene_attrs/gene_names")
    h5createDataset(hdf5_out, "raw_counts", dims=c(length(h5read(hdf5_in, "gene_attrs/gene_ids")),length(barcodes_subset)), chunk=c(100,100))
	h5write(h5read(hdf5_in, "raw_counts", index = list(NULL,cells_index)), file = hdf5_out, "raw_counts")
    
	h5closeAll()
}

create_hdf5_subset(barcodes_subset, opt$input_data, opt$output_subset)
