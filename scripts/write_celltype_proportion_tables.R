library(optparse)
library(SingleCellExperiment)
library(rhdf5)
library(dplyr)
library(knitr)
library(kableExtra)

option_list = list(
  make_option("--data_full", type = "character", help = "Path to RDS file with full CDS data object."),
  make_option("--data_subset", type = "character", help = "Path to RDS file with subsetted CDS data object."),
  make_option("--data_raw", type = "character", help = "Path to barcodes file from cellranger count output."),
  make_option("--data_subset_raw", type = "character", help = "Path to HDF5 file containing subsetted raw counts data."),
  make_option("--data_dict", type = "character", help = "Path to the ground-truth dictionary file mapping cell types barcodes."),
  make_option("--output", type = "character", help = "Directory to where final table(s) will be written."),
  make_option("--sample_name", type = "character", help = "Name of sample to prepend to output table file name.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

data_full_path <- opt$data_full
data_subset_path <- opt$data_subset
data_raw_path <- opt$data_raw
data_subset_raw_path <- opt$data_subset_raw
data_dictionary_path <- opt$data_dict

print(opt$data_full)
if(length(opt$data_full)){data_full_sce <- readRDS(data_full_path)}else{data_full_sce <- ""}
#data_raw_barcodes <- h5read(opt$data_raw,"cell_attrs/cell_names")
#data_raw_barcodes <- read.csv(opt$data_raw)[,1] # Use 10X barcodes
data_raw_barcodes <- tryCatch({
data_raw_barcodes <- h5read(opt$data_raw,"cell_attrs/cell_names") # Use rawCounts hdf5 file
}, error = function(e){
        data_raw_barcodes <- read.csv(opt$data_raw)[,1]
}
)

data_subset_sce <- readRDS(data_subset_path)

data_subset_raw_barcodes <- tryCatch({
data_subset_raw_barcodes <- h5read(opt$data_subset_raw,"cell_attrs/cell_names") # Use rawCounts hdf5 file
}, error = function(e){
	data_subset_raw_barcodes <- read.csv(opt$data_subset_raw)[,1]
}
)
data_dictionary <- read.table(data_dictionary_path, sep = '\t', header = TRUE)
print(data_raw_barcodes[1]);print(data_dictionary$Barcode[1])
if(nchar(as.character(data_raw_barcodes[1])) < nchar(as.character(data_dictionary$Barcode[1]))){
	print("Barcodes in data and in ground truth are of inequal lengths. Adjusting barcode lengths...")
	print(data_dictionary$Barcode[1])
	data_dictionary$Barcode <- substr(as.vector(data_dictionary$Barcode), 1 , nchar(data_raw_barcodes[1]))
	print(data_dictionary$Barcode[1])
}

barcode_match <- match(as.vector(data_raw_barcodes), as.vector(data_dictionary$Barcode))
#print(barcode_match)
total_raw_barcodes <- length(data_dictionary$Barcode)
total_data_barcodes <- length(data_raw_barcodes)
print(anyNA(barcode_match))
print(paste0("Barcodes matched: ", sum(lengths(na.omit(barcode_match)))))
print(paste0("% matched: ", sum(lengths(na.omit(barcode_match)))/ length(barcode_match) * 100))
print(total_raw_barcodes)
print(total_data_barcodes)

#data_raw_barcodes <- data.frame(as.vector(data_raw_barcodes))
#data_raw_barcodes$true_labels <- data_dictionary$Original_type[barcode_match]
data_raw_barcodes <- data_dictionary[barcode_match, ]

barcode_match <- match(as.vector(data_subset_raw_barcodes), as.vector(data_dictionary$Barcode))
#data_subset_raw_barcodes <- data.frame(as.vector(data_subset_raw_barcodes))
#data_subset_raw_barcodes$true_labels <- as.vector(data_dictionary$Original_type)[barcode_match]
data_subset_raw_barcodes <- data_dictionary[barcode_match, ]

print(head(data_raw_barcodes))
print(head(data_subset_raw_barcodes))

print(head(data_raw_barcodes$Original_type))
print(names(data_raw_barcodes))
print(length(data_raw_barcodes$Barcode))
print(length(data_subset_raw_barcodes$Barcode))

#sample_list <- list(
#  "_full_processed" = data_full_sce@metadata,
#  "_full_raw" = data_raw_barcodes,
#  "_subset_processed" = data_subset_sce@metadata,
#  "_subset_raw" = data_subset_raw_barcodes,
#  "_metadata" = data_dictionary
#)

#if(!length(opt$data_full)){sample_list <- sample_list[-1]}

if(!length(opt$data_full)){
    sample_list <- list(
  "_full_raw" = data_raw_barcodes,
  "_subset_processed" = data_subset_sce@metadata,
  "_subset_raw" = data_subset_raw_barcodes,
  "_metadata" = data_dictionary
)} else{
	    sample_list <- list(
  "_full_processed" = data_full_sce@metadata,
  "_full_raw" = data_raw_barcodes,
  "_subset_processed" = data_subset_sce@metadata,
  "_subset_raw" = data_subset_raw_barcodes,
  "_metadata" = data_dictionary
)}

######## LATEST VER
final_tables <- c()
tables <- list(
  #"ground_truth_major"="table(sample$ground_truth_major)",
  #"ground_truth_minor"="table(sample$ground_truth_minor)",
  "true_labels"="table(sample$true_labels)"
  )

dictionary_tables <- list(
  #"coerced_major_type"="table(sample$Coerced_major_type)",
  #"coerced_minor_type"="table(sample$Coerced_minor_type)",
  "original_labels" = "table(sample$Original_type)"
	)

for(sample_table in names(tables)){
  final_table <- data.frame()
  
  for(sample_name in names(sample_list)){
  print(sample_name)
  print(sample_table)
 	  if(grepl("raw", sample_name) || grepl("metadata", sample_name)){
		  print("Not HDF5/SCE file.")
    	sample <- data.frame(sample_list[sample_name], check.names = FALSE)
	print(names(sample))
	names(sample) <- c("Barcode", "Original_type", "Coerced_major_type", "Coerced_minor_type")
    	sample_label <- case_when(
  	 			sample_table == "ground_truth_major" ~ "coerced_major_type",
 		 			sample_table == "ground_truth_minor" ~ "coerced_minor_type",
 		 			sample_table == "true_labels" ~ "original_labels",
 		 			TRUE ~ as.character(sample_table)
			)
	freq_table <- eval(parse(text = dictionary_tables[sample_label]))
  	} else{
		print("HDF5/SCE file detected.")
    sample <- sample_list[sample_name]$"_"
    names(sample)
    freq_table <- eval(parse(text = tables[sample_table]))
    }
  	prop_table <- as.data.frame(round(prop.table(freq_table)*100, 1))
  	full_table <- data.frame(freq_table)
  	colnames(full_table) <- c("Cell_type","Frequency")
  	full_table$"Proportion" <- prop_table[,2]
  	full_table[,2] <- sapply(full_table$Frequency, as.numeric)
  	full_table[,3] <- sapply(full_table$Proportion, as.double)
 	 levels(full_table[,1]) <- c(levels(full_table[,1]), "Total")
	 totals <- c("Total"=c("Total",sum(full_table[,2]), sum(full_table[,3])))
	 print(totals)
 	 full_table[nrow(full_table)+1, ] <- totals
 	 rownames(full_table) <- c(rownames(full_table)[1:nrow(full_table)-1], "Total")
         names(full_table)[2:3] <- sapply(names(full_table)[2:3], paste0, sample_name)
	 print(full_table)

  	if(length(final_table) > 0){


		sym_diff <- function(a,b) setdiff(union(a,b), intersect(a,b))
		missing_cell_types <- sym_diff(final_table$Cell_type, full_table$Cell_type)
		for(cell_type in missing_cell_types){
			if(!cell_type %in% full_table$Cell_type)
			{ print(cell_type); print("Cell type not in current table.")
			#full_table[nrow(full_table)+1,] <- c(cell_type,rep(0, ncol(full_table)-1))
			all_cell_types <- c(as.vector(full_table$Cell_type),cell_type)
			full_table[nrow(full_table)+1,] <- c(cell_type,rep(0, ncol(full_table)-1))
			full_table$Cell_type <- as.factor(all_cell_types)
			}

		if(!cell_type %in% final_table$Cell_type)
                        { print(cell_type); print("Cell type not in previous table(s).")
                        #full_table[nrow(full_table)+1,] <- c(cell_type,rep(0, ncol(full_table)-1))
                        #final_table$Cell_type <- as.factor(c(as.vector(final_table$Cell_type), cell_type))
                        final_table[nrow(final_table)+1,] <- c(cell_type,rep(0, ncol(final_table)-1))
                        }
		}
		print("Merging...")
  		final_table <- merge(final_table, full_table)
    	} else {final_table <- rbind(final_table, full_table)}
     
   # print(final_table)
   
    #write.table(final_table, paste0(opt$output, opt$sample_name, "_", sample_table, "_table.txt"), quote = F, row.names = F, sep = "\t")
    #print("Table saved.")
  }

    colnames(final_table) <- sapply(strsplit(colnames(final_table), "_"), function(x) {
        tmp <- x[1]
        x[1] <- x[length(x)]
        x[length(x)] <- tmp
        x <- tolower(x)
        paste0(strsplit(x, "_"), collapse="_")
    })
    colnames(final_table)[1] <- "Cell_type"

    write.table(final_table, paste0(opt$output, opt$sample_name, "_", sample_table, "_table.txt"), quote = F, row.names = F, sep = "\t")
    print("Table saved. [1/3]")

    #metadata_index <- grep("metadata", sapply(strsplit(colnames(final_table), "_"), function(i) i[2]))
    #raw_index <- grep("raw", sapply(strsplit(colnames(final_table), "_"), function(i) i[3]))
    metadata_index <- grep("metadata", colnames(final_table))
    raw_index <- grep("raw", colnames(final_table))
    processed_index <- grep("processed", colnames(final_table))
    frequency_index <- grep("frequency", colnames(final_table))
    proportion_index <- grep("proportion", colnames(final_table))

    print(final_table)
    write.table(final_table[, c(1, metadata_index, raw_index, processed_index)], paste0(opt$output, opt$sample_name, "_", sample_table, "_table_v2.txt"), quote = F, row.names = F, sep = "\t")
    print("Table saved. [2/3]")

   write.table(final_table[, c(1, frequency_index, proportion_index)], paste0(opt$output, opt$sample_name, "_", sample_table, "_table_v3.txt"), quote = F, row.names = F, sep = "\t")
    print("Table saved. [3/3]") 
  	final_tables <- c(final_tables, list(final_table))
  }



#####
# FINAL OUTPUT
# 
#options("browser"="firefox")
print("Writing PDFs...")
for(i in seq(length(final_tables))){
  metadata_count <- sum(lengths(regmatches(colnames(final_tables[[i]]), gregexpr("metadata", colnames(final_tables[[i]])))))
  full_count <- sum(lengths(regmatches(colnames(final_tables[[i]]), gregexpr("full", colnames(final_tables[[i]])))))
  subset_count <- sum(lengths(regmatches(colnames(final_tables[[i]]), gregexpr("subset", colnames(final_tables[[i]])))))
  columns <- c(" "=1, "Full" = full_count, "Subset" = subset_count, "Metadata" = metadata_count)
  final_tables[[i]] %>%
    knitr::kable(caption = paste(opt$sample_name, "cell type proportions")) %>% 
      kable_styling(latex_options = "striped") %>%
	add_header_above(columns) %>%
          row_spec(nrow(final_tables[[i]]), bold = TRUE) %>%
  	    save_kable(paste0(opt$sample_name, "_", names(tables)[i], "_table.pdf"))
    	    print("Saved.")
  }
