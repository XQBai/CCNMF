
PLot_CNV_subclone_heatmap <-function(CNVsample, label, Gene_Order, index_remain_cells){
  
  Label <- matrix(0, nrow=dim(CNVsample)[2])
  Label[index_remain_cells] <- label
  Label <- Label + 1
  
  Gene_Order[['Gene_ID']] = as.character(Gene_Order[["gene"]])
  Gene_Order[['Gene_Chr']] = as.character(Gene_Order[["chr"]])
  row.names(Gene_Order) = Gene_Order[['Gene_ID']]
  
  CNVsample[["Gene_ID"]] = rownames(CNVsample)
  #there are genes in gistic output not represented in Gene_Order
  #and vice versa
  #length(CNVsample[["Gene_ID"]])
  Gene_Order = Gene_Order[Gene_Order[['Gene_ID']] %in% CNVsample[['Gene_ID']],]
  CNVsample = CNVsample[CNVsample[["Gene_ID"]] %in% Gene_Order[['Gene_ID']],]
  
  row_genes <- CNVsample[['Gene_ID']]
  CNVsample$Gene_ID <- NULL
  row.names(CNVsample) <- row_genes  

  ###bug
  CNVsample  <- CNVsample[Gene_Order[['Gene_ID']], ] #reorder Charlie_Gistic2_CNR 
  # CNVsample_noise <- CNVsample[, setdiff(CNV_barcode$V1, colnames(CNVsample1))] #subsect noisy cells 
  
  # CNVsample <- CNVsample[, colnames(CNVsample1)] #subset remain cells
  
  rownames(CNVsample) <- Gene_Order[['Gene_ID']]
  # rownames(CNVsample_noise) <- Gene_Order[['Gene_ID']]
  
  #mat = t(CNVsample)[,1:300] #reduce size for debugging
  mat = t(CNVsample)
  mat <- SortCells(t(mat), Label)
  mat[mat > 6] = 6
  
  Labels <- rep(paste0('s', sort(unique(Label)-1)), table(Label))
  Labels <- data.frame(cluster = Labels)
  
  gene_chr = rle(Gene_Order[["Gene_Chr"]][match(colnames(mat), Gene_Order[["Gene_ID"]])])
  gene_chr_mark = c(0, cumsum(gene_chr$lengths))[seq_along(gene_chr$lengths)]+1
  names(gene_chr_mark) = gene_chr$values
  Plot_CNV_heatmap(mat, Labels, gene_chr_mark, 'CNV_subclones.pdf')
}

SortCells <- function(tmp_scdna_matrix_arm, S1){
  ## Re-order cells according to the cluster labels
  clusterlabel <- c()
  barcodes <- c()
  
  tmp_matrix <- as.data.frame(matrix(0, dim(tmp_scdna_matrix_arm)[1], dim(tmp_scdna_matrix_arm)[2]))
  for (j in 1:max(S1)){
    clusterlabel <- c(clusterlabel, length(which(S1 == j)))
    barcodes <- c(barcodes, colnames(tmp_scdna_matrix_arm)[which(S1 == j)])
    if(j == 1){
      tmp_matrix[, 1:clusterlabel] = tmp_scdna_matrix_arm[, which(S1 == j)]
    }else{
      a = sum(clusterlabel[1:j-1]) + 1
      b = sum(clusterlabel[1:j])
      tmp_matrix[, a:b]=tmp_scdna_matrix_arm[, which(S1 == j)]
    }
  }
  tmp_matrix <- as.matrix(tmp_matrix)
  rownames(tmp_matrix) <- rownames(tmp_scdna_matrix_arm)
  colnames(tmp_matrix) <- barcodes
  tmp_matrix <- t(tmp_matrix)
  return(tmp_matrix)
}


Plot_CNV_heatmap <- function(mat, labels, gene_chr_mark, filename){
  
  col_fun = colorRamp2(c(0,1,2,3,4,5,6), c('blue', 'light blue', 'white', 'light salmon', 'coral', 'red', 'dark red'))
  ht_font = 40
  grid_height = 1
  
  if(length(unique(labels$cluster)) == 2){
    Label_color <- c("#00BFC4", "#F8766D")
    names(Label_color) <- sort(unique(labels$cluster))
  }else if(length(unique(labels$cluster)) >=3 & length(unique(labels$cluster)) <= 9){
    Label_color <- brewer.pal(length(unique(labels$cluster)), "Set1")
    names(Label_color) <- sort(unique(labels$cluster))
  }else if(length(unique(labels$cluster)) > 9){
    print("The number of clusters exceed 9, please add additional colorlabel for clusters.")
  }
  
  column_ha = ComplexHeatmap::columnAnnotation(Chr = ComplexHeatmap::anno_mark(at=gene_chr_mark[1:22],
                                                                               side="bottom",
                                                                               labels=names(gene_chr_mark)[1:22],
                                                                               labels_gp = grid::gpar(fontsize = ht_font)),
                                               gap=unit(1, 'mm'))
  
  
  row_ha = ComplexHeatmap::rowAnnotation(Label = labels$cluster,
                                         col = list(Label = Label_color),
                                         annotation_name_gp = grid::gpar(fontsize = ht_font),
                                         annotation_legend_param = list(title_gp = gpar(fontsize = ht_font, fontface = "bold"),
                                                                                   direction = "horizontal",
                                                                                    grid_height = grid::unit(grid_height,"inch"),
                                                                                    labels_gp = gpar(fontsize = ht_font)),
                                         annotation_name_side = list(Label = "bottom"))
  
  ht = ComplexHeatmap::Heatmap(mat,
                               name = "CNV", 
                               show_row_dend = FALSE,
                               show_column_dend = FALSE,
                               show_row_names = FALSE,
                               show_column_names = FALSE,
                               col = col_fun,
                               row_title_gp = gpar(col = c("blue", "red"), fontsize = ht_font),
                               bottom_annotation = column_ha,
                               column_names_centered = TRUE,
                               right_annotation = row_ha,
                               cluster_columns = FALSE,
                               cluster_rows = FALSE,
                               # column_km = 22,
                               # column_gap = unit(0, "in"),
                               border = TRUE,
                               heatmap_legend_param = list(title_gp = gpar(fontsize = ht_font, fontface = "bold"),
                                                           direction = "horizontal",
                                                           grid_height = grid::unit( grid_height,"inch" ),
                                                           labels_gp = gpar(fontsize = 0.8 * ht_font))
                                )
  
  pdf(file=filename, width=40, height=20)
  ComplexHeatmap::draw(ht, merge_legend = T, heatmap_legend_side = "bottom", annotation_legend_side = 'bottom')
  dev.off()
}




