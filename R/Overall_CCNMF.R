#' @description whole progress of running CCNMF
#' @import NMF
#' @import stats
#' @param CNVmatrix_input input cnv matrix
#' @param RNAmatrix_input inout rna matrix
#' @param ncluster the number of given clusters
#' @param CoupledMatrix the user-defined coupled matrix
#' @param initial_parameters the set of lambda1 and lambda2
#' @param lambda1 the weight of NMF
#' @param lambda2 the weight of coupled item
#' @param alpha hyper-parameter alpha
#' @param beta hyper-parameter beta
#' @param initial_coupling set up coupling matrix
#' @return list of results
#'
#' @export
Integrate_CCNMF <- function(CNVmatrix_input, RNAmatrix_input, ncluster,
                            CoupledMatrix = NULL,
                            initial_parameters = c('hyper-parameter', 'same-order', 'user-defined'),
                            lambda1 = 1, lambda2 = 2, alpha = NULL, beta = NULL,
                            initial_coupling = c('default', 'user-defined')){

  #utils::globalVariables(c('ncluster'))
  # Remove all 0 rows in CNVmatrix and RNAmatrix
  label_mean1 <- apply(CNVmatrix_input, 1, mean)
  index_mean1 <- which(label_mean1 != 0)
  CNVmatrix_input <- CNVmatrix_input[index_mean1, ]
  RNAmatrix_input <- RNAmatrix_input[index_mean1, ]
  label_mean2 <- apply(RNAmatrix_input, 1, mean)
  index_mean2 <- which(label_mean2 != 0)
  CNVmatrix_input <- CNVmatrix_input[index_mean2, ]
  RNAmatrix_input <- RNAmatrix_input[index_mean2, ]

  # Construct the identity matrix as copuled matrix, or use an informative matrix
  if(is.null(CoupledMatrix)){
    CoupledMatrix <- as.matrix(diag(dim(CNVmatrix_input)[1]))
  }else{
    CoupledMatrix <- CoupledMatrix
  }

  # Initialize lambda1 and lambda2 by seting the items in object function in the same order
  if(initial_parameters == 'hyper-parameter'){
    # Run NMF on scDNA matrix and scRNA matrix respectively, cost much time for large matrix
    resCNV <- NMF::nmf(CNVmatrix_input, ncluster)
    W10 <- resCNV@fit@W
    H10 <- resCNV@fit@H

    resRNA <- NMF::nmf(RNAmatrix_input, ncluster)
    W20 <- resRNA@fit@W
    H20 <- resRNA@fit@H
    par <- default_parameters_hyper(CNVmatrix_input, W10, H10, RNAmatrix_input, W20, H20,
                                    CoupledMatrix, alpha, beta)
    lambda1 <- par[1]
    lambda2 <- par[2]
  }else if(initial_parameters == 'same-order'){
    resCNV <- NMF::nmf(CNVmatrix_input, ncluster)
    W10 <- resCNV@fit@W
    H10 <- resCNV@fit@H

    resRNA <- NMF::nmf(RNAmatrix_input, ncluster)
    W20 <- resRNA@fit@W
    H20 <- resRNA@fit@H
    par <- default_parameters_order(CNVmatrix_input, W10, H10, RNAmatrix_input, W20, H20, CoupledMatrix)
    lambda1 <- par[1]
    lambda2 <- par[2]
  }else{
    lambda1 <- lambda1
    lambda2 <- lambda2
  }

  print(paste0('lambda1:', lambda1))
  print(paste0('lambda2:', lambda2))

  ResultsCCNMF <- run_CCNMF(ncluster, RNAmatrix_input, CNVmatrix_input, CoupledMatrix, lambda1, lambda2)

  H1 <- ResultsCCNMF[[1]]
  H2 <- ResultsCCNMF[[2]]

  ## Normalize the H1 and H2
  m1 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
  m2 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
  for (i in 1:ncluster){
    m1[i, i] <- sqrt(sum(H1[, i]^2))
    m2[i, i] <- sqrt(sum(H2[, i]^2))
  }

  H1 <- solve(m1) %*% H1
  H2 <- solve(m2) %*% H2
  S1 <- apply(H1, 2, function(x){which(x == max(x))})
  S2 <- apply(H2, 2, function(x){which(x == max(x))})

  ResultsCCNMF[[5]] <- S1
  ResultsCCNMF[[6]] <- S2

  # saveRDS(S1, file='S1.rds')
  # saveRDS(S2, file='S2.rds')
  return(ResultsCCNMF)
}

