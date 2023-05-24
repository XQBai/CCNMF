
#' @description Initialize the parameters lambda1 and lambda2 by hyper-parameter alpha and beta
#' @param CNVmatrix the copy number matrix
#' @param W10 W matrix of CNV's NMF
#' @param H10 H matrix of CNV's NMF
#' @param RNAmatrix gene expression matrix
#' @param W20 W matrix of RNA NMF
#' @param H20 H matrix of RNA's NMF
#' @param CoupledMatrix coupled Matrix
#' @param beta hyper-parmater
#' @param alpha hyper-parmater
#'
#' @export
default_parameters_hyper <- function(CNVmatrix, W10, H10, RNAmatrix, W20, H20, CoupledMatrix, beta, alpha){

  eps <- exp(-10)
  A <- CoupledMatrix
  if (is.null(beta) & is.null(alpha)){
    alpha = 0.001
    beta = 0.5}
  else if (is.numeric(beta) & is.null(alpha)){
    alpha = 0.001
  }

  r1 <- mean((apply(CNVmatrix %*% t(H10), MARGIN = 2, mean))/(apply(t(A) %*% W20, MARGIN=2, mean)))

  lambda2 <- 2*alpha*beta*r1

  r2 <- mean(apply(RNAmatrix %*% t(H20), MARGIN = 2, mean)/(apply(A %*% W10, MARGIN=2, mean)))

  lambda1 <- beta/(1-beta+eps)*r1/r2

  par <- c(lambda1, lambda2)
  return(par)
}

#' @description Initialize the parameters lambda1 and lambda2 by seting the items in object function in the same order
#' @param CNVmatrix the copy number matrix
#' @param W10 W matrix of CNV's NMF
#' @param H10 H matrix of CNV's NMF
#' @param RNAmatrix gene expression matrix
#' @param W20 W matrix of RNA NMF
#' @param H20 H matrix of RNA's NMF
#' @param CoupledMatrix coupled Matrix
#'
#' @export
default_parameters_order <- function(CNVmatrix,W10, H10, RNAmatrix, W20, H20, CoupledMatrix){

  E <- RNAmatrix
  O <- CNVmatrix
  A <- CoupledMatrix

  lambda1 <- norm((E - W20 %*% H20), type <- 'F')**2 / norm((O-W10 %*% H10), type <- 'F')**2
  ncluster <- dim(W10)[2]

  # Normalize the W10 and W20
  m1 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
  m2 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
  for (i in 1:ncluster){
    m1[i, i] <- sqrt(sum(W10[, i]^2))
    m2[i, i] <- sqrt(sum(W20[, i]^2))
  }

  W10 <- W10 %*% solve(m1)
  W20 <- W20 %*% solve(m2)

  lambda2 <-  norm((O-W10 %*% H10), type <- 'F')**2 / sum(diag(t(W20) %*% A %*% W10))

  lambda1 <- 10^(floor(log10(lambda1)))
  lambda2 <- 10^(floor(log10(lambda2))-3)
  par <- c(lambda1, lambda2)
  return(par)
}
