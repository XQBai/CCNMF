
#' Output the all visualization of results from CCNMF
#' Especially, the paired heatmap of differential genes and dimension reduction for both scRNA-seq and scDNA-seq data
#' @import ggplot2
#' @import cowplot
#' @import ggplotify
#'
#' @description Output the all visualization of results from CCNMF
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
#' @import ggplot2
#' @import grDevices
#'
#' @param CNVmatrix copy number matrix
#' @param RNAmatrix gene expression matrix
#' @param Result_CCNMF The result of CCNMF
#' @param DElist The list of differentiate genes
#'
#' @return The integrated figure
#' @export
PlotMainResult <- function(CNVmatrix, RNAmatrix, Result_CCNMF, DElist){

  if(is.list(DElist)){
    P <- DElist[[1]]
  }else if (is.character(DElist)){
    P <- DElist
  }
  DNAheat <- Plot_heatmap(CNVmatrix, Result_CCNMF[[5]], P, Datatype = 'DNA', title = 'Signature gene heatmap of scDNA-seq data')
  DNAheat <- as.ggplot(DNAheat)
  RNAheat <- Plot_heatmap(RNAmatrix, Result_CCNMF[[6]], P, Datatype = 'RNA', title = 'Signature gene heatmap of scRNA-seq data')
  RNAheat <- as.ggplot(RNAheat)

  H1 <- Result_CCNMF[[1]]
  H2 <- Result_CCNMF[[2]]
  S1 <- Result_CCNMF[[5]]
  S2 <- Result_CCNMF[[6]]
  DNAdim <- Plottsne(H1, S1, 'Tsne plot of scDNA-seq data', Datatype = 'scDNA-seq', need_PCA = FALSE)
  RNAdim <- Plottsne(H2, S2, 'Tsne plot of scRNA-seq data', Datatype = 'scRNA-seq', need_PCA = FALSE)

  myplot <- plot_grid(DNAdim, RNAdim, DNAheat, RNAheat, labels = c('A', 'B', 'C', 'D'), label_size = 12, scale = c(1, 1, 0.95, 0.95))
  dev.off()
  ggsave(filename ='allfigure.pdf', plot = myplot, width = 8.5, height = 6)
}


#' Plot the dimensional figure for scRNA-seq and scDNA-seq datasets
#' The rows in Data matrix represent the cells, the columns represent the genes/chr bins
#'
#' Plot the dimensional reduction figure by PCA
#' @import stats
#' @param Data the orgial matrix or H matrix when the raw matrix after NMF
#' @param label the clusters label of this figure which is the result of CCNMF
#' @param Datatype The type of input Data matrix, 'scRNA-seq' or 'scDNA-seq'
#' @param title the title of the figure
#' @param need_PCA logic
#'
#' @return The pdf figure based on PCA
#' @export
Plotpca <- function(Data, label, title, Datatype = 'scRNA-seq'){
  pca <- prcomp(t(Data))
  myplot <- Plot(pca$x[, 1:2], as.character(label), title, 'pca 1', 'pca 2')
  return(myplot)
}

#' Plot the dimensional reduction figure by tsne
#' Before tsne, the PCA is applied to the original Data, then select the top 15 PCs and apply
#' tsne to these components.
#' @import stats
#' @param Data the orgial matrix or H matrix when the raw matrix after NMF
#' @param label the clusters label of this figure which is the result of CCNMF
#' @param Datatype The type of input Data matrix, 'scRNA-seq' or 'scDNA-seq'
#' @param title the title of the figure
#' @param need_PCA logic
#'
#' @return The pdf figure based on tsne
#' @export
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
#' @import stats
#' @param Data the orgial matrix or H matrix when the raw matrix after NMF
#' @param label the clusters label of this figure which is the result of CCNMF
#' @param Datatype The type of input Data matrix, 'scRNA-seq' or 'scDNA-seq'
#' @param need_PCA Logic
#'
#' @return The pdf figure based on umap
#' @export
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
#' @import ggplot2
#' @param Data the input data which need to be plotted, which can be raw matrix or H matrix concluded by NMF
#' @param Cluster the clusters label
#' @param title the title of the figure, which is a string
#' @param labelx X-coordinate of name, such as 'pca 1', 'tsne 1', 'umap 1'
#' @param labely Y-coordinate of name, such as 'pca 2', 'tsne 2', 'umap 2'
#'
#' @return A pdf file
#' @export
Plot <- function(Data, Cluster, title, labelx,  labely){
  Data <- as.data.frame(Data)
  colnames(Data) <- c('V1', 'V2')
  ncluster <- max(Cluster)
  Clones <- paste0('C', Cluster, sep='')
  myplot <- ggplot(data = Data, ggplot2::aes(x = V1, y=V2, color = Clones)) +
    geom_point(size = 1) +
    labs(title= title, x = labelx, y= labely, fill = 'Clones') +
    scale_fill_discrete(name='Clones') +
    theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 11,  color = 'black', face = 'bold', vjust=0.5, hjust = 0.5)) +
    theme(axis.title.y = element_text(size = 11, color = 'black', face = 'bold', vjust=0.5, hjust = 0.5)) +
    theme(axis.text.x = element_text(size = 7, color = 'black', face = 'bold')) +
    theme(axis.text.y = element_text(size = 7, color = 'black', face = 'bold')) +
    theme(legend.text= element_text(size=7,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
    theme(legend.title = element_text(size=7,color="black", face= "bold"))
  if(length(unique(Cluster)) == 2){
    Label_color <- c('C1' = "#00BFC4", 'C2' = "#F8766D")
  }else if(length(unique(Cluster)) >=3 & length(unique(Cluster)) <= 9){
    Label_color <- brewer.pal(length(unique(Cluster)), "Set1")
    label_cluster <- paste0('C', Cluster)
    names(Label_color) <- sort(unique(label_cluster))
  }else if(length(unique(Cluster)) > 9){
    print("The number of clusters exceed 9, please add additional colorlabel for clusters.")
  }

  myplot <- myplot + scale_colour_manual(values = Label_color)
  return(myplot)
}

#' Plot the heatmap for differential genes in scRNA-seq data
#' @import RColorBrewer
#' @import grDevices
#' @import pheatmap
#'
#' @param Data scRNA-seq gene expression matrix
#' @param label the clustering label of scRNA-seq data
#' @param Datatype scDNAseq or scRNA-seq
#' @param P the p_values matrix by t test the DE genes among clusters in scRNA-seq data
#' @param title a string, default as 'The heatmap of differential expression in scRNA-seq data
#'
#' @return D_assign: a matrix that rows are the clusters and cloumns are the DE genes of all clusters
#' @return A pdf file which is the heatmap figure
#' @export
Plot_heatmap <- function(Data, label, P, Datatype = 'DNA', title = 'The heatmap of differential expression in scRNA-seq data'){

  # if input P is a dataframe, which means the p-values for each cluster
  # otherwise, P is the vector of differential genes, whose format is character.
  if(is.data.frame(P)){
    index <- apply(P, 2, function(x){order(x, decreasing=F)[1:10]})
    DEgene <- c()
    for (i in 1: dim(index)[2]){
      DEgene <- c(DEgene, rownames(P)[index[, i]])
    }
  }else if(is.character(P)){
    DEgene <- P
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

  annotation_row = data.frame(Clones = factor(rep(paste('C', 1:max(label), sep = ''), cluster_num)))
  rownames(annotation_row) = cell_name

  if(Datatype == 'DNA'){
    myplot <- pheatmap::pheatmap(D_assign, color=colorRampPalette(rev(c("red","white","blue")))(102), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = F, fontsize = 7, annotation_row = annotation_row, angle_col = "45",annotation_legend = FALSE, main = title)
  }else{
    myplot <- pheatmap::pheatmap(D_assign, color=colorRampPalette(rev(c("red","white","blue")))(102), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = F, fontsize = 7, annotation_row = annotation_row, angle_col = "45",annotation_legend = FALSE, main = title)
  }
  return(myplot)
}


