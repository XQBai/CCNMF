#' Output the all visualization of results from CCNMF
#' Especially, the paired heatmap of differential genes and dimension reduction for both scRNA-seq and scDNA-seq data
#' @description Output the all visualization of results from CCNMF
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
#' @importFrom ggplot2 ggsave
#' @importFrom cowplot plot_grid
#' @param CNVmatrix copy number matrix
#' @param RNAmatrix gene expression matrix
#' @param Result_CCNMF The result of CCNMF
#'
#' @return The integrated figure
#' @export
PlotMainResult <- function(CNVmatrix, RNAmatrix, Result_CCNMF){

  H1 <- Result_CCNMF[[1]]
  H2 <- Result_CCNMF[[2]]
  S1 <- paste0('C', Result_CCNMF[[5]])
  S2 <- paste0('R', Result_CCNMF[[6]])
  DNAdim <- Plottsne(H1, S1, 'tSNE plot of scDNA-seq data', Datatype = 'scDNA-seq', need_PCA = FALSE)
  RNAdim <- Plottsne(H2, S2, 'tSNE plot of scRNA-seq data', Datatype = 'scRNA-seq', need_PCA = FALSE)

  myplot <- plot_grid(DNAdim, RNAdim, labels = c('A', 'B'), label_size = 12, scale = c(1, 1))
  ggsave(filename ='CCNMF_Tsne_plot.pdf', plot = myplot, width = 10, height = 4.5)
}


#' Plot the dimensional figure for scRNA-seq and scDNA-seq datasets
#' The rows in Data matrix represent the cells, the columns represent the genes/chr bins
#'
#' Plot the dimensional reduction figure by PCA
#' @importFrom stats prcomp
#' @param Data the orgial matrix or H matrix when the raw matrix after NMF
#' @param label the clusters label of this figure which is the result of CCNMF
#' @param Datatype The type of input Data matrix, 'scRNA-seq' or 'scDNA-seq'
#' @param title the title of the figure
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
#' @importFrom stats prcomp
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
      tsne <- Rtsne(pca$x[, 1:50], dims = 2, perplexity = 30, check_duplicates = FALSE, verbose = TRUE, max_iter = 500)
    }
  }else{
    tsne <- Rtsne(t(Data), dims = 2, perplexity = 30, check_duplicates = FALSE, verbose = TRUE, max_iter = 500)
  }
  myplot <- Plot(tsne$Y, as.character(label), title, 'Tsne_1', 'Tsne_2')
  return(myplot)
}

#' Plot the dimensional reduction figure by umap
#' Firstly, we apply PCA for the high-dimensional data, then use Umap to the top 15 PCs.
#' @importFrom stats prcomp
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
      umap <- umap(pca$x[, 1:100], n_neighbors = 15, learning_rate = 0.5, init = "random")
    }
  }else{
    umap <- umap(t(Data), n_neighbors = 15, learning_rate = 0.5, init = "random")
  }
  myplot <- Plot(umap, as.character(label), Datatype, 'UMAP 1', 'UMAP 2')
  return(myplot)
}

#' The function of plotting figure
#' @import ggplot2
#' @importFrom rlang .data
#' @param Data the input data which need to be plotted, which can be raw matrix or H matrix concluded by NMF
#' @param Clones the clusters label
#' @param title the title of the figure, which is a string
#' @param labelx X-coordinate of name, such as 'pca 1', 'tsne 1', 'umap 1'
#' @param labely Y-coordinate of name, such as 'pca 2', 'tsne 2', 'umap 2'
#'
#' @return A pdf file
#' @export
Plot <- function(Data, Clones, title, labelx,  labely){
  Data <- as.data.frame(Data)
  colnames(Data) <- c('V1', 'V2')
  # ncluster <- max(Cluster)
  # Clones <- paste0('C', Cluster, sep='')
  myplot <- ggplot(data = Data, ggplot2::aes(x = .data$V1, y= .data$V2, color = Clones)) +
    geom_point(size = 0.5) +
    theme_classic() +
    labs(title= title, x = labelx, y= labely, fill = 'Clones') +
    scale_fill_discrete(name='Clones') +
    theme(plot.title = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 15,  color = 'black', face = 'bold', vjust=0.5, hjust = 0.5)) +
    theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold', vjust=0.5, hjust = 0.5)) +
    theme(axis.text.x = element_text(size = 10, color = 'black', face = 'bold')) +
    theme(axis.text.y = element_text(size = 10, color = 'black', face = 'bold')) +
    theme(legend.text= element_text(size=10,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
    theme(legend.title = element_text(size=10,color="black", face= "bold"))
  if(length(unique(Clones)) == 2){
    # Label_color <- c('C1' = "#e41a1c", 'C2' = "#377eb8")
    Label_color <- c("#e41a1c", "#377eb8")
    names(Label_color) <- sort(unique(Clones))
  }else if(length(unique(Clones)) >=3 & length(unique(Clones)) <= 9){
    Label_color <- brewer.pal(length(unique(Clones)), "Set1")
    # label_cluster <- paste0('C', Cluster)
    names(Label_color) <- sort(unique(Clones))
  }else if(length(unique(Clones)) > 9){
    print("The number of clusters exceed 9, please add additional colorlabel for clusters.")
  }

  myplot <- myplot + scale_colour_manual(values = Label_color)
  return(myplot)
}




