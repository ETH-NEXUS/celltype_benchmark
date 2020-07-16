# Adapted from (Abdelaal et al. 2019) https://github.com/tabdelaal/scRNAseq_Benchmark - script(s): evaluate.R 

option_list = list(
  make_option("--sce", type = "character", help = "Path to SCE object w/ ground-truth annotations."),
  make_option("--predicted_labels_path", type = "character", help = "Path to csv file w/ list of predicted labels produced by classifier."),
  make_option("--output_dir", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sample_name", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--method_name", type = "character", help = "Cell typing method identifier. Attached to each output name.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


evaluate_classifier <- function(sce, predicted_labels_path, indices = NULL){
  sce <- readRDS(opt$sce)
  true_labels <- sce@meta.data$ground_truth_major
  rm(sce)
  #true_labels <- unlist(read.csv(true_labels_path))
  predicted_labels <- unlist(read.csv(predicted_labels_path))
  
  if (! is.null(indices)){
    true_labels <- true_labels[indices]
    predicted_labels <- predicted_labels[indices]
  }
  
  unique_true <- unlist(unique(true_labels))
  unique_predicted <- unlist(unique(predicted_labels))
  unique_all <- unique(c(unique_true,unique_predicted))
  
  confusion_matrix <- table(true_labels,predicted_labels)
  pop_size <- rowSums(confusion_matrix)
  predicted_labels = gsub('Node..','Node',predicted_labels)
  confusion_matrix_f1_score <- table(true_labels,predicted_labels,exclude = c('unassigned','Unassigned','Unknown','rand','Node','ambiguous','unknown'))

  f1_score <- vector()
  sum_acc <- 0
  
  for (i in c(1:length(unique_true))){
    find_label = colnames(confusion_matrix_f1_score) == row.names(confusion_matrix_f1_score)[i]
    if(sum(find_label)){
      prec <- confusion_matrix_f1_score[i,find_label] / colSums(confusion_matrix_f1_score)[find_label]
      rec <- confusion_matrix_f1_score[i,find_label] / rowSums(confusion_matrix_f1_score)[i]
      if (prec == 0 || rec == 0){
        f1_score[i] = 0
      } else{
        f1_score[i] <- (2*prec*rec) / (prec + rec)
      }
      sum_acc <- sum_acc + confusion_matrix_f1_score[i,find_label]
    } else {
      f1_score[i] = 0
    }
  }
  
  pop_size <- pop_size[pop_size > 0]
  names(f1_score) <- names(pop_size)
  median_f1_score <- median(f1_score)
  total <- length(predicted_labels)
  num_unlab <- sum(predicted_labels == 'unassigned') + sum(predicted_labels == 'Unassigned') + sum(predicted_labels == 'rand') + sum(predicted_labels == 'Unknown') + sum(predicted_labels == 'unknown') + sum(predicted_labels == 'Node') + sum(predicted_labels == 'ambiguous')
  per_unlab <- num_unlab / total
  acc <- sum_acc/sum(confusion_matrix_f1_score)
  
  result <- list(conf = confusion_matrix, median_f1_score = median_f1_score, f1_score = f1_score, accuracy = acc, percent_unlabeled = per_unlab, population_size = pop_size)
  
  return(result)
}

results <- evaluate_classifier(opt$sce, opt$predicted_labels_path)
write.csv(results$conf, file.path(opt$output_dir,opt$sample_name, ".confusion_", paste0(opt$method_name, ".csv")))
write.csv(results$f1_score, file.path(opt$output_dir,opt$sample_name, ".f1-score_", paste0(opt$method_name, ".csv")))
write.csv(results$population_size, file.path(opt$output_dir, opt$sample_name,".population-size_", paste0(opt$method_name, ".csv")))
df <- data.frame(results[c("median_f1_score", "accuracy", "percent_unlabeled")])
write.csv(df, file.path(opt$output_dir, opt$sample_name, ".summary_", paste0(opt$method_name, ".csv")))
