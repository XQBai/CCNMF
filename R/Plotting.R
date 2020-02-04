
#' Output the all visualization of results from CCNMF
#' Especially, the paired heatmap of differential genes and dimension reduction for both scRNA-seq and scDNA-seq data
PlotMainResult <- function(CNVmatrix, RNAmatrix, Result_CCNMF, DElist){

  #ncluster <- max(Result_CCNMF[[5]])
  DNAheat <- Plot_heatmap(CNVmatrix, Result_CCNMF[[5]], DElist[[1]], Datatype = 'DNA', title = 'The heatmap of scDNA-seq data')
  DNAheat <- as.ggplot(DNAheat)
  RNAheat <- Plot_heatmap(RNAmatrix, Result_CCNMF[[6]], DElist[[1]], Datatype = 'RNA', title = 'The heatmap of scRNA-seq data')
  RNAheat <- as.ggplot(RNAheat)
  #DNAheat <- DNAheat + scale_colour_manual(values = hue_pal()(ncluster))
  #RNAheat <- DNAheat + scale_colour_manual(values = hue_pal()(ncluster))
  H1 <- Result_CCNMF[[1]]
  H2 <- Result_CCNMF[[2]]
  S1 <- Result_CCNMF[[5]]
  S2 <- Result_CCNMF[[6]]
  DNAdim <- Plottsne(H1, S1, 'Tsne for scDNA-seq data', Datatype = 'scDNA-seq', need_PCA = FALSE)
  RNAdim <- Plottsne(H2, S2, 'Tsne for scRNA-seq data', Datatype = 'scRNA-seq', need_PCA = FALSE)

  #DNAdim <- DNAdim + scale_colour_manual(values = hue_pal()(ncluster))
  #RNAdim <- RNAdim + scale_colour_manual(values = hue_pal()(ncluster))

  myplot <- plot_grid(DNAheat, RNAheat, DNAdim, RNAdim, labels = c('A', 'B', 'C', 'D'), label_size = 12, scale = c(0.95, 0.95, 1, 1))
<<<<<<< HEAD
  dev.off()
=======
>>>>>>> c773e3c756507e2ce1cd3d284df1d6b8f44dbb4b
  ggsave(file='allfigure.pdf', plot = myplot, width = 8.5, height = 6)
}


#' Plot the dimensional figure for scRNA-seq and scDNA-seq datasets
#' The rows in Data matrix represent the cells, the columns represent the genes/chr bins
#'
#' Plot the dimensional reduction figure by PCA
#'
#' @param Data the orgial matrix or H matrix when the raw matrix after NMF
#' @param label the clusters label of this figure which is the result of CCNMF
#' @param Datatype The type of input Data matrix, 'scRNA-seq' or 'scDNA-seq'
#'
#' @return The pdf figure based on PCA
Plotpca <- function(Data, label, title, Datatype = 'scRNA-seq'){
  pca <- prcomp(t(Data))
  myplot <- Plot(pca$x[, 1:2], as.character(label), title, 'pca 1', 'pca 2')
  return(myplot)
}

#' Plot the dimensional reduction figure by tsne
#' Before tsne, the PCA is applied to the original Data, then select the top 15 PCs and apply
#' tsne to these components.
#'
#' @param Data the orgial matrix or H matrix when the raw matrix after NMF
#' @param label the clusters label of this figure which is the result of CCNMF
#' @param Datatype The type of input Data matrix, 'scRNA-seq' or 'scDNA-seq'
#'
#' @return The pdf figure based on tsne
Plottsne <- function(Data, label, title, Datatype = 'scRNA-seq', need_PCA = TRUE){
  if (need_PCA == TRUE){
    pca <- prcomp(t(Data))
    if(dim(pca$x)[2] < 15){
      tsne <- Rtsne(pca$x, dims = 2, perplexity = 30, check_duplicates = FALSE, verbose = TRUE, max_iter = 500)
    }else{
      tsne <- Rtsne(pca$x[, 1:15], dims = 2, perplexity = 30, check_duplicates = FALSE, verbose = TRUE, max_iter = 500)
    }
  }else{
    tsne <- Rtsne(t(Data), dims = 2, perplexity = 30, check_duplicates = FALSE, verbose = TRUE, max_iter = 500)
  }
  myplot <- Plot(tsne$Y, as.character(label), title, 'tsne 1', 'tsne 2')
  return(myplot)
}

#' Plot the dimensional reduction figure by umap
#' Firstly, we apply PCA for the high-dimensional data, then use Umap to the top 15 PCs.
#'
#' @param Data the orgial matrix or H matrix when the raw matrix after NMF
#' @param label the clusters label of this figure which is the result of CCNMF
#' @param Datatype The type of input Data matrix, 'scRNA-seq' or 'scDNA-seq'
#'
#' @return The pdf figure based on umap
Plotumap <- function(Data, label, Datatype = 'scRNA-seq', need_PCA = TRUE){
  if (need_PCA == TRUE){
    pca <- prcomp(t(Data))
    if(dim(pca$x)[2] < 15){
      umap <- umap(pca$x, n_neighbors = 15, learning_rate = 0.5, init = "random")
    }else{
      umap <- umap(pca$x[, 1:15], n_neighbors = 15, learning_rate = 0.5, init = "random")
    }
  }else{
    umap <- umap(t(Data), n_neighbors = 15, learning_rate = 0.5, init = "random")
  }
  myplot <- Plot(umap, as.character(label), paste0('Umap for ', Datatype), 'umap 1', 'umap 2')
  return(myplot)
}

#' The function of plotting figure
#' @param Data the input data which need to be plotted, which can be raw matrix or H matrix concluded by NMF
#' @param label the clusters label
#' @param title the title of the figure, which is a string
#' @param labelx X-coordinate of name, such as 'pca 1', 'tsne 1', 'umap 1'
#' @param labely Y-coordinate of name, such as 'pca 2', 'tsne 2', 'umap 2'
#' @param savefig the name of the pdf file, such as "figure.pdf"
#'
#' @return A pdf file
Plot <- function(Data, Cluster, title, labelx,  labely){
  Data <- as.data.frame(Data)
  colnames(Data) <- c('V1', 'V2')
  ncluster <- max(Cluster)
  Cluster <- paste0('C', Cluster, sep='')
  myplot <- ggplot(data = Data, aes(x = V1, y=V2, color = Cluster)) +
    geom_point(size = 1) +
    labs(title= title, x = labelx, y= labely, fill = 'Cluster') +
    scale_fill_discrete(name='Clone') +
    theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 11,  color = 'black', face = 'bold', vjust=0.5, hjust = 0.5)) +
    theme(axis.title.y = element_text(size = 11, color = 'black', face = 'bold', vjust=0.5, hjust = 0.5)) +
    theme(axis.text.x = element_text(size = 7, color = 'black', face = 'bold')) +
    theme(axis.text.y = element_text(size = 7, color = 'black', face = 'bold')) +
    theme(legend.text= element_text(size=7,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
    theme(legend.title = element_text(size=7,color="black", face= "bold"))

  if (ncluster == 2){
    myplot <- myplot + scale_colour_manual(values = c('C1' = "#00BFC4", 'C2' = "#F8766D"))
  }else if (ncluster == 3){
    myplot <- myplot + scale_colour_manual(values = c('C1' = "#00BA38", 'C2' = "#F8766D", 'C3' = "#619CFF"))
  }else if (ncluster == 4){
    myplot <- myplot + scale_colour_manual(values = c('C1' = "#A3A500", 'C2' = "#F8766D", 'C3' = "#00BF7D", 'C4' = "#E76BF3"))
  }
  return(myplot)
}


#' Plot the heatmap for differential genes in scRNA-seq data
#'
#' @param Data scRNA-seq gene expression matrix
#' @param label the clustering label of scRNA-seq data
#' @param P the p_values matrix by t test the DE genes among clusters in scRNA-seq data
#' @param title a string, default as 'The heatmap of differential expression in scRNA-seq data
#'
#' @return D_assign: a matrix that rows are the clusters and cloumns are the DE genes of all clusters
#' @return A pdf file which is the heatmap figure
Plot_heatmap <- function(Data, label, P, Datatype = 'DNA', title = 'The heatmap of differential expression in scRNA-seq data'){

  index <- apply(P, 2, function(x){order(x, decreasing=F)[1:7]})
  DEgene <- c()
  for (i in 1: dim(index)[2]){
    DEgene <- c(DEgene, rownames(P)[index[, i]])
  }
  D <- Data[unique(DEgene), ]
  cluster_num <- c()
  cell_name <- c()
  D_assign = as.data.frame(matrix(0, dim(D)[1], dim(D)[2]))
  for (j in 1:max(label)){
    cluster_num <- c(cluster_num, length(which(label == j)))
    cell_name <-c(cell_name, colnames(D)[which(label ==j)])
    if (j == 1){
      D_assign[, 1: cluster_num[1]] =  D[, which(label == j)]
    }
    else{
      a = sum(cluster_num[1:j-1]) + 1
      b = sum(cluster_num[1:j])
      D_assign[, a : b] = D[, which(label == j)]
    }
  }
  rownames(D_assign) <- rownames(D)
  colnames(D_assign) <- cell_name
  D_assign <- as.data.frame(t(D_assign))

  annotation_row = data.frame(Cluster = factor(rep(paste('C', 1:max(label), sep = ''), cluster_num)))
  rownames(annotation_row) = cell_name

  if(Datatype == 'DNA'){
    myplot <- pheatmap(D_assign, color=colorRampPalette(rev(c("red","white","blue")))(102), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = F, fontsize = 7, annotation_row = annotation_row, angle_col = "45",annotation_legend = FALSE, main = title)
  }else{
    myplot <- pheatmap(D_assign, color=colorRampPalette(rev(c("red","white","blue")))(102), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = F, fontsize = 7, annotation_row = annotation_row, angle_col = "45",annotation_legend = FALSE, main = title)
  }
  return(myplot)
}


