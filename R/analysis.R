
#' @description Analysis the results of CCNMF
#' @param RNAmatrix gene expression matrix
#' @param CNVmatrix copy number matrix
#' @param ncluster the number of clusters
#' @param lambda1 paramaters
#' @param lambda2 paramaters
#' @param QC Logic
#' @param GeneFilter Logic
#' @param reference_name name of reference
#'
#' @export
AnalyzeResult <- function(RNAmatrix, CNVmatrix, ncluster = 3, lambda1, lambda2, QC = TRUE, GeneFilter=TRUE, reference_name = 'hg38'){

  if (QC == TRUE){
    RNAmatrix <- QualityControl(RNAmatrix)
  }
  if (GeneFilter == TRUE){
    RNAmatrix <- gene_filter(RNAmatrix)
  }
  Result <- Estimate_Coupled_matrix(RNAmatrix, CNVmatrix, reference_name = reference_name)
  E <- Result[[2]]
  P <- Result[[1]]
  A <- Result[[3]]
  ResultsCCNMF <- run_CCNMF(ncluster, E, P, A, lambda1, lambda2)
  return(list(Result, ResultsCCNMF))
}


#' @description Find the paired dosage effect genes
#' @import minerva
#' @param CNVmatrix copy number matrix
#' @param RNAmatrix gene expression matrix
#' @param CNVDE differenitial expressed genes
#' @param RNADE differential expressed genes
#' @param S1 clustering labels of DNA
#' @param S2 clustering labels of RNA
#' @param index logic
#'
#' @export
DosageQualify <- function(CNVmatrix, RNAmatrix, CNVDE, RNADE, S1, S2, index=FALSE){
  ncluster <- max(S1)
  Dosage <- list()
  MICmatrix <- matrix(0, nrow = dim(RNAmatrix)[1], ncol = ncluster)
  for (i in 1:ncluster){
    L1 <- apply(CNVmatrix, 1, function(x){x[which(S1 == i)]})
    L2 <- apply(RNAmatrix, 1, function(x){x[which(S2 == i)]})

    if (dim(L1)[1] >= dim(L2)[1]){
      L2 <- apply(L2, 2, function(x){sample(x, dim(L1)[1], replace = TRUE)})
    }else{
      L1 <- apply(L1, 2, function(x){sample(x, dim(L2)[1], replace = TRUE)})
    }

    for (j in 1:dim(RNAmatrix)[1]){
      MICmatrix[j , i] <- mine(L1[, j], L2[, j])$MIC
    }

    RNA_P <- RNADE[[1]]
    number <- 100
    RNAtop100 <- rownames(RNAmatrix)[order(RNA_P[ , i])[1: number]]
    MICtop100 <- rownames(RNAmatrix)[order(MICmatrix[ , i])[1: number]]
    Dosagene <- intersect(MICtop100, RNAtop100)
    if (index == TRUE){
      CNV_P <- CNVDE[[1]]
      CNVtop100 <- rownames(CNVmatrix)[order(CNV_P[ , i])[1: number]]
      Dosagene <- intersect(MICtop100, intersect(CNVtop100, RNAtop100))
    }
    Dosage[[i]] <- Dosagene
  }
  # Dosage <- as.data.frame(Dosage)
  return(Dosage)
}

#' @import minerva
#' @param CNVmatrix copy number matrix
#' @param RNAmatrix gene expression matrix
#' @param S1 clustering labels of DNA
#' @param S2 clustering labels of RNA
#' @export
ComputeMIC <- function(CNVmatrix, RNAmatrix, S1, S2){
  ncluster <- max(S1)
  Dosage <- list()
  MICmatrix <- matrix(0, nrow = dim(RNAmatrix)[1], ncol = ncluster)
  for (i in 1:ncluster){
    L1 <- apply(CNVmatrix, 1, function(x){x[which(S1 == i)]})
    L2 <- apply(RNAmatrix, 1, function(x){x[which(S2 == i)]})

    if (dim(L1)[1] >= dim(L2)[1]){
      L2 <- apply(L2, 2, function(x){sample(x, dim(L1)[1], replace = TRUE)})
    }else{
      L1 <- apply(L1, 2, function(x){sample(x, dim(L2)[1], replace = TRUE)})
    }

    for (j in 1:dim(RNAmatrix)[1]){
      MICmatrix[j , i] <- mine(L1[, j], L2[, j])$MIC
    }
  }
  return(MICmatrix)
}

#' Differential expression for different clusters in specific clusters
#' Selct the differential gene among clusters by t test
#' @import stats
#' @param Data the scRNA-seq matrix
#' @param label the clusters label of the imput data
#' @param i the i-th cluster
#'
#' @return P the p_values matrix, the rows correspond to the clusters, the cloumns correspond to the genes
#' @return DEgenes, the differential genes for each clusters. the each row is the DE genes for each cluster.
#' @export
DiffExpCluster <- function(Data, label, i){

  P <- matrix(0, dim(Data)[1], max(label))
  Cluster_specific_genes <- c()
  for (j in 1:max(label)){
    ###sapply
    if(j != i){
      p_val <- unlist(
        x = sapply(X = 1:nrow(Data),
          FUN = function(x){
            t.test(x = Data[x, which(label == i)], y = Data[x, which(label== j)])$p.value
          }
        )
      )
      P[, j] <- p_val
    }
  }
  P <- P[, -i]
  P <- as.data.frame(P, row.names = rownames(Data), sep='')
  #genesnames <- apply(P, 2, function(x){rownames(Data)})
  for (k in 1:dim(P)[2]){
    Cgenes <- rownames(Data)[order(P[, k])[1:10]]
    Cluster_specific_genes <- c(Cluster_specific_genes, Cgenes)
  }
  #DEgenes <- apply(P, 2, function(x){rownames(Data)[which(x < 0.05)]})
  return(list(P, Cluster_specific_genes))
}

#' Differential expression for different clusters in scRNA-seq data
#' Selct the differential gene among clusters by t test
#' @import stats
#' @param Data the scRNA-seq matrix
#' @param label the clusters label of the imput data
#'
#' @return P the p_values matrix, the rows correspond to the clusters, the cloumns correspond to the genes
#' @return DEgenes, the differential genes for each clusters. the each row is the DE genes for each cluster.
#' @export
DiffExp <- function(Data, label){

  std <- apply(Data, 1, sd)
  Data <- Data[which(std != 0), ]
  P <- matrix(0, dim(Data)[1], max(label))
  for (j in 1:max(label)){
    ###sapply
    p_val <- unlist(
      x = sapply(X = 1:nrow(Data),
        FUN = function(x){
          t.test(x = Data[x, which(label == j)], y = Data[x, which(label!= j)])$p.value
        }
      )
    )
    P[, j] <- p_val
  }
  P <- as.data.frame(P, row.names = rownames(Data), sep='')
  genesnames <- apply(P, 2, function(x){rownames(Data)})
  DEgenes <- apply(P, 2, function(x){rownames(Data)[which(x < 0.05)]})
  return(list(P, genesnames, DEgenes))
}

#' @description Find high variable genes for each cluster by vst
#' @param Data the scRNA-seq matrix
#' @param label the clusters label of the imput data
#' @export
FindHighVar <- function(Data, label, number){
  HighVar <- matrix(0, number, max(label))
  for (j in 1:max(label)){
    RNACj <- Data[, which(label == j)]
    m <- apply(RNACj, 1, mean)
    v <- apply(RNACj, 1, var)
    vst <- v/m
    i <- order(vst, decreasing = T)[1:number]
    HighVar[ , j] <- i
  }
  return(HighVar)
}

