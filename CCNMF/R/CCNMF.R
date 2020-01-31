
#' Run CCNMF function
#' Run \code{CCNMF} accross a range pf initializations
#'
#' @param ncluster the number of subpopulations/clones seted by users
#' @param RNAmatrix A matrix of gene expression in scRNA-seq
#' @param CNVmatrix A matrix of copy number in scDNA-seq
#' @param CoupledMatrix A coupled matrix between genes and cnv bins
#' @param lambda1 A hyper-parameter which is seted by users
#' @param lambda2 A hyper-parameter which is seted by users
#'
#' @details
#' This function is core of CCNMF and the parameters lambda1 and lambda2 control the objective function, meanwhile, the two hyper-parameters can be
#' initially assigend.
#'
#' @export
#'
#' @return The list(H1, H2, S1, S2)
#'
#' @examples
#' library(NMF)
#' RNAmatrix <-
#' CNVmatrix <-
#' CoupledMatrix <-
#' Result <- run_CCNMF(ncluster, RNAmatrix, CNVmatrix, CoupledMatrix, lambda1, lambda2)

run_CCNMF <- function(ncluster, RNAmatrix, CNVmatrix, CoupledMatrix, lambda1, lambda2, mu=1){

  start <- proc.time()
  print('Initializing NMF for CNV matrirx...')
  resCNV <- nmf(CNVmatrix, ncluster)
  W10 <- resCNV@fit@W
  H10 <- resCNV@fit@H

  err2 <- numeric(50)
  print('Initializing NMF for scRNA matrirx...')
  # for (i in 1:50){
  #   resRNA <- nmf(RNAmatrix, ncluster, seed <- i)
  #   W20 <- resRNA@fit@W
  #   H20 <- resRNA@fit@H
  #   err2[i] <- norm((RNAmatrix -W20 %*% H20 ), type <- 'F')
  # }
  resRNA <- nmf(RNAmatrix, ncluster)
  W20 <- resRNA@fit@W
  H20 <- resRNA@fit@H

  print('Initializing the parameters lambda1, lambda2 and mu...')

  print('Start Coupled NMF')

  E <- RNAmatrix
  O <- CNVmatrix
  A <- CoupledMatrix
  W1 <- W10
  W2 <- W20
  H1 <- H10
  H2 <- H20

  #mu <- 1
  eps <- 0.0001

  print('Iterating coupledNMF...')
  maxiter <- 500
  err <- 1
  terms <- numeric(maxiter)
  it <- 1
  ## The objective function ==========================================================================
  terms[it] <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2) + lambda1/2 * (norm(E - W2 %*% H2)**2)
  - lambda2 * sum(diag(t(W2) %*% A %*% W1)) + mu * (norm(W1, type <- 'F') ** 2 + norm(W2,type <- 'F') ** 2)
  while (it <<- maxiter & err > 1e-6){
    it <- it + 1
    T1 <- as.matrix(lambda2 * (t(A) %*% W2))
    T1[which(T1 < 0)] <- 0
    W1 <- W1 * ( O %*% t(H1) + T1) / (eps + W1 %*% H1 %*% t(H1) + 2 * mu * W1)
    H1 <- H1 * (t(W1) %*% O / (t(W1) %*% W1 %*% H1))

    T2 <- as.matrix((lambda2 / lambda1) * (t(A) %*% W1))
    T2[which(T2 < 0)] <- 0
    W2 <- W2 * ( E %*% t(H2) + T2) / (eps + W2 %*% H2 %*% t(H2) + 2*mu/lambda1 * W2)
    H2 <- H2 * (t(W2) %*% E)/(t(W2) %*% W2 %*% H2)

    ## Uniformize H1 and H2 =============================================================================
    m1 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
    m2 <- matrix(0, nrow <- ncluster, ncol <- ncluster)

    # for (i in 1:ncluster){
    #   m1[i, i] <- sqrt(sum(H1[i,]^2))
    #   m2[i, i] <- sqrt(sum(H2[i, ]^2))
    # }
    #
    # W1 <- W1 %*% m1
    # W2 <- W2 %*% m2
    # H1 <- solve(m1) %*% H1
    # H2 <- solve(m2) %*% H2

    terms[it] <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2) + lambda1/2 * (norm(E - W2 %*% H2)**2)
    - lambda2 * sum(diag(t(W2) %*% A %*% W1)) + mu * (norm(W1, type <- 'F') ** 2 + norm(W2, type <- 'F') ** 2)
    err <- abs(terms[it] - terms[it - 1]) / abs(terms[it-1] + eps)

    eq1 <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2)
    eq2 <- lambda1/2 * (norm(E - W2 %*% H2)**2)
    eq3 <- lambda2 * sum(diag(t(W2) %*% A %*% W1))
    eq4 <- mu * (norm(W1, type <- 'F') ** 2 + norm(W2, type <- 'F') ** 2)
    eq5 <- terms[it]

    # W1 <- W1 %*% solve(m1)
    # W2 <- W2 %*% solve(m2)
}

for (i in 1:ncluster){
  m1[i, i] <- sqrt(sum(H1[i,]^2))
  m2[i, i] <- sqrt(sum(H2[i, ]^2))
}
  # W1 <- W1 %*% m1
  # W2 <- W2 %*% m2
  H1 <- solve(m1) %*% H1
  H2 <- solve(m2) %*% H2
  S1 <- apply(H1, 2, function(x){which(x == max(x))})
  S2 <- apply(H2, 2, function(x){which(x == max(x))})
  end <- proc.time()
#
  H1 <- m1 %*% H1
  H2 <- m2 %*% H2

  print(paste0('eq1 ', round(eq1)))
  print(paste0('eq2 ', round(eq2)))
  print(paste0('eq3 ', round(eq3)))
  print(paste0('eq4 ', round(eq4)))
  print(paste0('eq5 ', round(eq5)))
  value <- rep(1, 5)
  value[1] <- round(eq1)
  value[2] <- round(eq2)
  value[3] <- round(eq3)
  value[4] <- round(eq4)
  value[5] <- round(eq5)

  print(paste0('Run time:', (end-start)[3][[1]], 'seconds'))
  return(list(H1, H2, W1, W2, S1, S2, value))
}





