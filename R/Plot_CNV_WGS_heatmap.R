#' @import RColorBrewer
#' @param CNVsample copy number data matrix with all genes
#' @param label subclones labeling of CCNMF
#' @param Gene_order dataframe of genes locs
#' @param index_remain_cells cell index of non-replicating cells
#' @return the figure of subclones heatmap
#' @export
PLot_CNV_subclone_heatmap <-function(CNVsample, label, Gene_Order, index_remain_cells){

  Label <- matrix(0, nrow=dim(CNVsample)[2])
  Label[index_remain_cells] <- label
  Label <- Label + 1

  Gene_Order[['Gene_ID']] = as.character(Gene_Order[["gene"]])
  Gene_Order[['Gene_Chr']] = as.character(Gene_Order[["chr"]])
  row.names(Gene_Order) = Gene_Order[['Gene_ID']]

  CNVsample <- as.data.frame(CNVsample)
  CNVsample[["Gene_ID"]] = rownames(CNVsample)

  #there are genes in gistic output not represented in Gene_Order
  #and vice versa
  Gene_Order = Gene_Order[Gene_Order[['Gene_ID']] %in% CNVsample[['Gene_ID']],]
  CNVsample = CNVsample[CNVsample[["Gene_ID"]] %in% Gene_Order[['Gene_ID']],]

  row_genes <- CNVsample[['Gene_ID']]
  CNVsample$Gene_ID <- NULL
  row.names(CNVsample) <- row_genes

  ###bug
  CNVsample  <- CNVsample[Gene_Order[['Gene_ID']], ] #reorder Charlie_Gistic2_CNR

  rownames(CNVsample) <- Gene_Order[['Gene_ID']]

  mat = CNVsample
  # Recorder cells in each subclones based on the distance to normal CNV profile
  normal_CNV <- matrix(2, nrow = 1, ncol = dim(CNVsample)[1])
  d <- as.matrix(dist(t(cbind(t(as.matrix(normal_CNV)), mat))))[1, ]
  d <- d[1:dim(mat)[2]+1]
  
  mat <- SortCells(as.matrix(mat), Label, d)
  mat[mat > 6] = 6

  Labels <- rep(paste0('C', sort(unique(Label)-1)), table(Label))
  Labels <- data.frame(cluster = Labels)

  gene_chr = rle(Gene_Order[["Gene_Chr"]][match(colnames(mat), Gene_Order[["Gene_ID"]])])
  gene_chr_mark = c(0, cumsum(gene_chr$lengths))[seq_along(gene_chr$lengths)]+1
  names(gene_chr_mark) = gene_chr$values

  Plot_CNV_heatmap(mat, Labels, gene_chr_mark, 'Subclones_no_noise.pdf')
}

#' @param mat the matrix of segments-bins
#' @param labels subclones index for cells
#' @param hc results of hierarchical clustering
#' @return list of ordered matrix and cluster label
#'
Sort_subclones <- function(mat, labels, hc){
  order <- paste0('C', hc$order)
  index <- c()
  for (i in order){
    index_sub <- which(labels$cluster == i)
    index <- c(index, index_sub)
  }
  labels$cluster <- labels$cluster[index]
  mat <- mat[index, ]
  return(list(mat, labels))
}

#' @param tmp_scdna_matrix the matrix of cells
#' @param Corr the distance to normal cell
#' @return Ordered cells matrix
#' 
SortCells_in_subclone <- function(tmp_scdna_matrix, Corr){
  
  tmp_scdna_matrix <- as.data.frame(t(tmp_scdna_matrix))
  tmp_scdna_matrix$corr <- Corr
  tmp_scdna_matrix <- tmp_scdna_matrix %>% arrange(corr)
  tmp_scdna_matrix$corr <- NULL
  tmp_scdna_matrix <- t(tmp_scdna_matrix)
  return(tmp_scdna_matrix)
}

#' @param tmp_scdna_matrix_arm scDNA matrix
#' @param S1 labels of subclonses
#' @param Corr the distance to normal cell
#' @return Ordered cells based on per subclones
#' @export
SortCells <- function(tmp_scdna_matrix_arm, S1, Corr){
  ## Re-order cells according to the cluster labels
  clusterlabel <- c()
  barcodes <- c()
  
  tmp_matrix <- as.data.frame(matrix(0, dim(tmp_scdna_matrix_arm)[1], dim(tmp_scdna_matrix_arm)[2]))
  for (j in 1:length(unique(S1))){
    clusterlabel <- c(clusterlabel, length(which(S1 == j)))
    barcodes <- c(barcodes, colnames(tmp_scdna_matrix_arm)[which(S1 == j)])
    if(j == 1){
      tmp_matrix[, 1:clusterlabel] = tmp_scdna_matrix_arm[, which(S1 == j)]
      tmp_matrix[, 1:clusterlabel] <- SortCells_in_subclone(tmp_matrix[, 1:clusterlabel], Corr[which(S1 == j)])
    }else{
      a = sum(clusterlabel[1:j-1]) + 1
      b = sum(clusterlabel[1:j])
      tmp_matrix[, a:b]=tmp_scdna_matrix_arm[, which(S1 == j)]
      tmp_matrix[, a:b] <- SortCells_in_subclone(tmp_matrix[, a:b], Corr[which(S1 == j)])
    }
  }
  tmp_matrix <- as.matrix(tmp_matrix)
  rownames(tmp_matrix) <- rownames(tmp_scdna_matrix_arm)
  colnames(tmp_matrix) <- barcodes
  tmp_matrix <- t(tmp_matrix)
  return(tmp_matrix)
}


#' @import RColorBrewer
#' @import ComplexHeatmap
#' @import circlize
#' @import ggplot2
#' @import grDevices
#' @import grid
#' @param mat matrix of copy number
#' @param labels subclones clustering labels
#' @param gene_chr_mark chromosome regions
#' @param filename the name of output figure
#' @return the pdf figure
#' @export
Plot_CNV_heatmap <- function(mat, labels, gene_chr_mark, filename){

  ht_font = 40
  grid_height = 1
  col_fun = colorRamp2(c(0,1,2,3,4,5,6), c('blue', 'light blue', 'white', 'light salmon', 'coral', 'red', 'dark red'))
  
  
  #Label_color <- c('C1' = "#00BA38", 'C2' = "#F8766D", 'C3' = "#619CFF")
  #Label_color <- c('C1' = "#00BFC4", 'C2' = "#F8766D")
if(length(unique(labels$cluster)) > 9){
  mycol = brewer.pal(9, 'Set1')
  Label_color <- c(mycol, brewer.pal(12, 'Set3')[10:12], "#666666")
  names(Label_color) <- sort(unique(labels$cluster))
}else if(length(unique(labels$cluster)) >= 3 & length(unique(labels$cluster)) <= 9){
  Label_color <- brewer.pal(length(unique(labels$cluster)), "Set1")
  names(Label_color) <- sort(unique(labels$cluster))
}else if(length(unique(labels$cluster)) == 2){
  Label_color <- c("red", "green")
  names(Label_color) <- sort(unique(labels$cluster))
}else if(length(unique(labels$cluster)) == 1){
  Label_color <- c('red')
  names(Label_color) <- sort(unique(labels$cluster))
}

  column_ha = HeatmapAnnotation(Chr = ComplexHeatmap::anno_mark(at=gene_chr_mark[1:22],
                                                                side="bottom",
                                                                labels=names(gene_chr_mark)[1:22],
                                                                labels_gp = grid::gpar(fontsize = ht_font)))

  ht1 = Heatmap(mat,name = 'CNV', na_col = 'white', show_row_dend = FALSE,
                show_column_dend = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                col = col_fun,
                border = TRUE,
                bottom_annotation = column_ha,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                row_gap = grid::unit(0.3, "in"),
                height = grid::unit(18, 'in'),
                width = grid::unit(36, 'in'),
                heatmap_legend_param = list(title_gp = gpar(fontsize = ht_font, fontface = "bold"),
                                            direction = "horizontal",
                                            grid_width = grid::unit(grid_height, 'inch'),
                                            grid_height = grid::unit(grid_height,"inch" ),
                                            labels_gp = gpar(fontsize = 0.8 * ht_font)))

  ht2 = Heatmap(labels$cluster,name = 'Subclone', col = Label_color, show_row_dend = FALSE,
                show_column_dend = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                border = TRUE,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                row_gap = grid::unit(0.3, "in"),
                width = grid::unit(0.5, 'in'),
                heatmap_legend_param = list(title_gp = gpar(fontsize = ht_font, fontface = "bold"),
                                            direction = "horizontal",
                                            grid_width = grid::unit(grid_height, 'inch'),
                                            grid_height = grid::unit(grid_height,"inch" ),
                                            labels_gp = gpar(fontsize = 0.8 * ht_font)))
  ht_list = ht1 + ht2

  pdf(file=filename, width=40, height=20)
  ComplexHeatmap::draw(ht_list,
                       ht_gap = unit(0.5, 'cm'),
                       heatmap_legend_side = "left",
                       annotation_legend_side = 'left',
                       merge_legend = T
  )
  chr_line <- gene_chr_mark/dim(mat)[2]
  chr_line[1] <- 0
  chr_line <- c(chr_line, 1)
  lapply(chr_line,
         FUN = function(p){
           decorate_heatmap_body("CNV", {
             grid.lines(c(p, p), c(0, 1), gp = gpar(lty = 1, lwd = 3))
           }, slice = 1)
         }
  )
  dev.off()
}



