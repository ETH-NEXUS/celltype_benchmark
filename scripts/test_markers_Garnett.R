###########################
###########################
# Test defined marker usefuleness using Garnett
library(org.Hs.eg.db)
library(garnett)

data_path <- "/cluster/project/nexus/benchmarking/celltyping/test_pilot/cellranger_run/zheng_sorted_merged/"

pbmc_cds <- load_cellranger_data(data_path)

# Basic pre-processing that NEEDS to be done
process_cds_monocle <- function(cds){
  cds <- preprocess_cds(cds, num_dim = 100) #PCA by default
  cds <- reduce_dimension(cds) #UMAP by default
  return(cds)
}

pbmc_cds <- process_cds_monocle(pbmc_cds)

marker_file_path <- c(
  "/cluster/project/nexus/benchmarking/celltyping/garnett_files/marker_files/hsPBMC_markers.txt", 
  "/cluster/project/nexus/benchmarking/celltyping/marker_files/PBMC/PBMC_Garnett_markers.txt"
  )

marker_check <- check_markers(cds = pbmc_cds, marker_file = marker_file_path[2], 
                              db = org.Hs.eg.db, classifier_gene_id_type = "SYMBOL",
                              cds_gene_id_type = "ENSEMBL")
#plot_markers(marker_check)

marker_check_list <- split(marker_check, marker_check$cell_type)
plot_markers <- rbind(marker_check_list[1:12])

ok_markers <- marker_check[marker_check$summary=="Ok",]
median(as.matrix(na.omit(marker_check$marker_score)))
high_markers <- ok_markers[ok_markers$marker_score > median(as.matrix(na.omit(marker_check$marker_score))),]
top_markers <- ok_markers[ok_markers$marker_score > 1,]

marker_check_zheng <- check_markers(cds = pbmc_cds, 
                                    marker_file = "/cluster/project/nexus/benchmarking/celltyping/marker_files/PBMC/PBMC_Garnett_markers_Zheng.txt", 
                              db = org.Hs.eg.db, classifier_gene_id_type = "SYMBOL",
                              cds_gene_id_type = "ENSEMBL")

plot_markers(marker_check_zheng[marker_check_zheng$summary=="Ok",])
####################TESTING##############################
plot_markers_DEBUG <- function (marker_check_df, amb_marker_cutoff = 0.5, label_size = 2) 
{
  assertthat::assert_that(is.data.frame(marker_check_df))
  assertthat::assert_that(sum(!c("marker_gene", "cell_type", 
                                 "nominates", "total_nominated", "most_overlap", "ambiguity", 
                                 "marker_score", "summary") %in% names(marker_check_df)) == 
                            0, msg = paste("marker_check_df must be the output of", 
                                           "the check_markers function. Input is", "missing key columns"))
  labeldf <- data.frame(cell_type = marker_check_df$cell_type, 
                        cell_type_label = paste0(marker_check_df$cell_type, ": ", 
                                                 marker_check_df$total_nominated), total_nominated = marker_check_df$total_nominated)
  labeldf <- labeldf[order(labeldf$total_nominated), ]
  labeldf <- labeldf[!duplicated(labeldf$cell_type), ]
  labeldf$cell_type <- as.factor(labeldf$cell_type)
  marker_check_df$cell_type <- as.factor(marker_check_df$cell_type)
  marker_check_df$cell_type <- factor(marker_check_df$cell_type, 
                                      labels = as.character(labeldf[order(labeldf$cell_type), 
                                                                    ]$cell_type_label))
  marker_check_df$marker_gene <- as.factor(marker_check_df$marker_gene)
  marker_check_df$tempy <- as.factor(paste(marker_check_df$marker_gene, 
                                           marker_check_df$cell_type, sep = ">"))
  
  marker_check_df1 <- #12 types
  marker_check_df2 <- #12 types
  marker_check_df3 <- #13 types
  ggplot2::ggplot(marker_check_df, ggplot2::aes(x = ambiguity, 
                                                y = forcats::fct_reorder2(tempy, cell_type, -marker_score), 
                                                fill = 100 * nominates/total_nominated)) + ggplot2::geom_point(ggplot2::aes(size = marker_score), 
                                                                                                               color = "black", pch = 21, data = marker_check_df, stroke = 0.1) + 
    ggplot2::geom_text(ggplot2::aes(x = 0.4, label = ifelse(is.na(ambiguity), 
                                                            as.character(summary), "")), color = "firebrick4", 
                       size = label_size) + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(ambiguity > 
                                                                                                   amb_marker_cutoff, as.character(paste0("High overlap with\n", 
                                                                                                                                          most_overlap)), "")), color = "black", point.padding = 0.01, 
                                                                     size = label_size, min.segment.length = ggplot2::unit(0, 
                                                                                                                           "lines"), segment.size = 0.2) + ggplot2::facet_grid(cell_type ~ 
                                                                                                                                                                                 ., scales = "free", space = "free_y", labeller = ggplot2::label_value) + 
    ggplot2::scale_size(name = "Marker\nscore", range = c(1, 
                                                          3)) + viridis::scale_fill_viridis(position = "top", 
                                                                                            name = "% of\nassigned") + ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 0.5, 
                                                                                                                                                                       barheight = 4)) + ggplot2::xlim(c(0, 1)) + ggplot2::ylab("") + 
    ggplot2::xlab("Ambiguity") + ggplot2::annotate("segment", 
                                                   x = -Inf, xend = Inf, y = Inf, yend = Inf, color = "grey93") + 
    ggplot2::annotate("segment", x = -Inf, xend = Inf, y = -Inf, 
                      yend = -Inf, color = "grey93") + ggplot2::theme_classic() + 
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0), 
                   strip.background = ggplot2::element_rect(color = "grey93", 
                                                            fill = "grey93"), legend.key.height = ggplot2::unit(0.2, 
                                                                                                                "cm")) + ggplot2::scale_y_discrete(labels = function(x) {
                                                                                                                  stringr::str_split_fixed(x, ">", 2)[, 1]
                                                                                                                })
}