#' Run CCNMF function
#' Run \code{CCNMF} accross a range pf initializations
#' @importFrom preprocessCore normalize.quantiles.robust
#' @importFrom RcppHungarian HungarianSolver
#' @importFrom stats runif
#' @importFrom lsa cosine
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
#' @return The list(H1, H2, S1, S2)
#'
#' @export
run_CCNMF <- function(ncluster, RNAmatrix, CNVmatrix, CoupledMatrix, lambda1, lambda2){

  tolx <- 1e-5
  tolfun <- 1e-6
  rate <- 1

  start <- proc.time()
  print('Initializing W10 and H10...')
  ngenes <- dim(CNVmatrix)[1]
  ncells <- dim(CNVmatrix)[2]
  W10 <- matrix(stats::runif(ngenes*ncluster, 0, 1), ngenes, ncluster)
  H10 <- matrix(stats::runif(ncluster * ncells, 0, 1), ncluster, ncells)

  print('Initializing W20 and H20 matrirx...')
  ngenes <- dim(RNAmatrix)[1]
  ncells <- dim(RNAmatrix)[2]
  W20 <- matrix(stats::runif(ngenes*ncluster, 0, 1), ngenes, ncluster)
  H20 <- matrix(stats::runif(ncluster * ncells, 0, 1), ncluster, ncells)
  # Normalize the W10 and W20
  m1 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
  m2 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
  for (i in 1:ncluster){
    m1[i, i] <- sqrt(sum(W10[, i]^2))
    m2[i, i] <- sqrt(sum(W20[, i]^2))
  }

  W10 <- W10 %*% solve(m1)
  W20 <- W20 %*% solve(m2)

  print('Initializing the parameters lambda1, lambda2...')

  print('Start Coupled NMF')
  E <- RNAmatrix
  O <- CNVmatrix
  A <- CoupledMatrix

  eps <- 0.0001
  err <- 1
  print('Iterating coupledNMF...')
  maxiter <- 1000
  terms <- numeric(maxiter)
  it <- 1
  dnorm0 <- norm((O-W10 %*% H10), type <- 'F')**2 + lambda1 * (norm(E - W20 %*% H20)**2)

  tmp_it <- c(1)
  eq10 <- 0.5*(norm((O - W10 %*% H10), type <- 'F')**2)
  eq20 <- lambda1/2 * (norm(E - W20 %*% H20)**2)
  eq30 <- lambda2 * sum(diag(t(W20) %*% A %*% W10))
  score1<- eq10 + eq20 - eq30
  tmp_score <- score1

  ## The objective function ==========================================================================
  txtfile <- file('output_ADMM_assign.txt', 'w')
  # & err > 1e-6
  while (it <= maxiter){
    it <- it + 1
    nassign <- 0

    T1 <- as.matrix(lambda2 * (t(A) %*% W20))
    T1[which(T1 < 0)] <- 0
    numer <- t(W10) %*% O
    H1 <- H10 * (numer / (t(W10) %*% W10 %*% H10) + eps)
    H1[which(H1 <0)] <- 0
    numer <- O %*% t(H1)
    mu11 <- diag(apply(W10 * (numer + T1), MARGIN = 2, sum))
    mu12 <- diag(apply(W10 * (W10 %*% (H1 %*% t(H1))), MARGIN = 2, sum))
    W1 <- W10 * (numer + T1 + W10 %*% mu12) / (eps + W10 %*% H1 %*% t(H1) + W10 %*% mu11)
    W1[which(W1 <0)] <- W1

    T2 <- as.matrix((lambda2 / (lambda1 + eps)) * (A %*% W1))
    T2[which(T2 < 0)] <- 0
    numer <- t(W20) %*% E
    H2 <- H20 * (numer)/(t(W20) %*% W20 %*% H20 + eps)
    H2[which(H2 < 0 )] <- 0
    numer <- E%*%t(H2)
    mu21 <- diag(apply(W20*(numer + T2), MARGIN = 2, sum))
    mu22 <- diag(apply(W20*(W20 %*% (H2 %*% t(H2))), MARGIN = 2, sum))
    W2 <- W20 * (numer + T2 + W20 %*% mu22) / (eps + W20 %*% H2 %*% t(H2) + W20 %*% mu21)
    W2[which(W2 < 0)] <- 0

    dnorm1 <- norm((O-W1 %*% H1), type <- 'F')**2 + lambda1 * (norm(E - W2 %*% H2)**2)
    dw1 <- max(abs(W1-W10))/(max(abs(W10)) + eps)
    dh1 <- max(abs(H1-H10))/(max(abs(H10)) + eps)
    dw2 <- max(abs(W2-W20))/(max(abs(W20)) + eps)
    dh2 <- max(abs(H2-H20))/(max(abs(H20)) + eps)
    delta <- max(c(dw1, dh1, dw2, dh2))
    if (it>1){
      if(delta <= tolx){
        print(paste0('delta=', delta, " is smaller"))
        break
      }
      else if (abs(dnorm1 - dnorm0) < 0.01){
        print(paste0('Stop at the optimal solution'))
        break
        }
      else if (it == maxiter){
        break
      }
    }

    W10 <- W1
    W20 <- W2
    H10 <- H1
    H20 <- H2


    if (it %% 50 == 0 & nassign == 0){
      # print('start')
      S10 <- apply(H10, 2, function(x){which(x == max(x))})
      S20 <- apply(H20, 2, function(x){which(x == max(x))})
      FC1 <- matrix(0, nrow <- ngenes, ncol <- ncluster)
      FC2 <- matrix(0, nrow <- ngenes, ncol <- ncluster)
      # print('before')
      for (i in 1:ncluster){
        FC1[, i] <- apply(O, 1, function(x){sum(x[which(S10 == i)] > 2)}) / apply(O, 1, function(x){sum(x[which(S10 != i)] > 2)} + 1) * (sum(S10 == i)/sum(S10 != i) + 1)
        FC2[, i] <- apply(E, 1, function(x){sum(x[which(S20 == i)] > 0)}) / apply(E, 1, function(x){sum(x[which(S20 != i)] >0)} + 1) * (sum(S20 == i)/sum(S20 != i) + 1)
      }
      WP1 <- preprocessCore::normalize.quantiles.robust(FC1, use.median = TRUE)
      WP2 <- preprocessCore::normalize.quantiles.robust(FC2, use.median = TRUE)
      # print('after')
      S <- t(WP2) %*% A %*% WP1
      # print("run S")
      S <- 1/(1 + ((S-min(S)))/(max(S) - min(S)))
      err <- (dnorm0-dnorm1)/dnorm0

      print(paste0("Interation: ", it))
      # print(paste0('dnorm0 ',  round(dnorm0, 2)))
      # print(paste0('dnorm1 ', round(dnorm1, 2)))
      # print(paste0("error ", (dnorm0-dnorm1)/dnorm0))
      #
      nassign <- nassign + 1

      eq1 <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2)
      eq2 <- lambda1/2 * (norm(E - W2 %*% H2)**2)
      eq3 <- lambda2 * sum(diag(t(W2) %*% A %*% W1))
      score1<- eq1 + eq2 - eq3

      tmp_it <- c(tmp_it, it)
      tmp_score <- c(tmp_score, score1)

      print(paste0('objective function ', round(score1, 2)))

    }
    dnorm0 <- dnorm1
  }
  W1 <- W10
  H1 <- H10
  W2 <- W20
  H2 <- H20

  # Score <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2) + lambda1/2 * (norm(E - W2 %*% H2)**2) - lambda2 * sum(diag(t(W2) %*% A %*% W1))
  S10 <- apply(H1, 2, function(x){which(x == max(x))})
  S20 <- apply(H2, 2, function(x){which(x == max(x))})

  #### Reassignment
  FC1 <- matrix(0, nrow <- ngenes, ncol <- ncluster)
  FC2 <- matrix(0, nrow <- ngenes, ncol <- ncluster)
  for (i in 1:ncluster){
    FC1[, i] <- apply(O, 1, function(x){mean(x[which(S10 == i)] )})
    FC2[, i] <- apply(E, 1, function(x){mean(x[which(S20 == i)] )})
  }

  ### Calculate the cosine similarity
  similarity <- matrix(0, nrow = ncluster, ncol = ncluster)
  for (i in 1:ncluster){
    for (j in 1:ncluster){
      similarity[i, j] <- round(lsa::cosine(FC1[, i], FC2[, j]), 3)
    }
  }

  S = 1 - similarity
  assignment <- HungarianSolver(S)$pairs
  W2 <- W2[, assignment[, 2]]
  H2 <- H2[assignment[, 2], ]
  S1 <- apply(H1, 2, function(x){which(x == max(x))})
  S2 <- apply(H2, 2, function(x){which(x == max(x))})

  eq1 <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2)
  eq2 <- lambda1/2 * (norm(E - W2 %*% H2)**2)
  eq3 <- lambda2 * sum(diag(t(W2) %*% A %*% W1))
  score <- eq1 + eq2 - eq3
  print(paste0('eq1 ', round(eq1)))
  print(paste0('eq2 ', round(eq2)))
  print(paste0('eq3 ', round(eq3)))
  out <- paste("Iteration", c(as.character(it), as.character(eq1), as.character(eq2), as.character(eq3)))
  writeLines(out, txtfile)
  close(txtfile)
  end <- proc.time()
  print(paste0('Run time:', (end-start)[3][[1]], 'seconds'))

  write.csv(data.frame(iteration = tmp_it,
                       score = tmp_score), file = 'covengence_table.txt')

  return(list(H1, H2, W1, W2, S1, S2))
}

#' Run CCNMF function
#' Run \code{CCNMF} accross a range pf initializations
#' @importFrom preprocessCore normalize.quantiles.robust
#' @importFrom RcppHungarian HungarianSolver
#' @importFrom stats runif
#' @importFrom lsa cosine
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
#' @return The list(H1, H2, S1, S2)
#'
#' @export
run_CCNMF <- function(ncluster, RNAmatrix, CNVmatrix, CoupledMatrix, lambda1, lambda2){

  tolx <- 1e-5
  tolfun <- 1e-6
  rate <- 1

  start <- proc.time()
  print('Initializing W10 and H10...')
  ngenes <- dim(CNVmatrix)[1]
  ncells <- dim(CNVmatrix)[2]
  W10 <- matrix(stats::runif(ngenes*ncluster, 0, 1), ngenes, ncluster)
  H10 <- matrix(stats::runif(ncluster * ncells, 0, 1), ncluster, ncells)

  print('Initializing W20 and H20 matrirx...')
  ngenes <- dim(RNAmatrix)[1]
  ncells <- dim(RNAmatrix)[2]
  W20 <- matrix(stats::runif(ngenes*ncluster, 0, 1), ngenes, ncluster)
  H20 <- matrix(stats::runif(ncluster * ncells, 0, 1), ncluster, ncells)
  # Normalize the W10 and W20
  m1 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
  m2 <- matrix(0, nrow <- ncluster, ncol <- ncluster)
  for (i in 1:ncluster){
    m1[i, i] <- sqrt(sum(W10[, i]^2))
    m2[i, i] <- sqrt(sum(W20[, i]^2))
  }

  W10 <- W10 %*% solve(m1)
  W20 <- W20 %*% solve(m2)

  print('Initializing the parameters lambda1, lambda2...')

  print('Start Coupled NMF')
  E <- RNAmatrix
  O <- CNVmatrix
  A <- CoupledMatrix

  eps <- 0.0001
  err <- 1
  print('Iterating coupledNMF...')
  maxiter <- 1000
  terms <- numeric(maxiter)
  it <- 1
  dnorm0 <- norm((O-W10 %*% H10), type <- 'F')**2 + lambda1 * (norm(E - W20 %*% H20)**2)

  tmp_it <- c(1)
  eq10 <- 0.5*(norm((O - W10 %*% H10), type <- 'F')**2)
  eq20 <- lambda1/2 * (norm(E - W20 %*% H20)**2)
  eq30 <- lambda2 * sum(diag(t(W20) %*% A %*% W10))
  score1<- eq10 + eq20 - eq30
  tmp_score <- score1

  ## The objective function ==========================================================================
  txtfile <- file('output_ADMM_assign.txt', 'w')
  # & err > 1e-6
  while (it <= maxiter){
    it <- it + 1
    nassign <- 0

    T1 <- as.matrix(lambda2 * (t(A) %*% W20))
    T1[which(T1 < 0)] <- 0
    numer <- t(W10) %*% O
    H1 <- H10 * (numer / (t(W10) %*% W10 %*% H10) + eps)
    H1[which(H1 <0)] <- 0
    numer <- O %*% t(H1)
    mu11 <- diag(apply(W10 * (numer + T1), MARGIN = 2, sum))
    mu12 <- diag(apply(W10 * (W10 %*% (H1 %*% t(H1))), MARGIN = 2, sum))
    W1 <- W10 * (numer + T1 + W10 %*% mu12) / (eps + W10 %*% H1 %*% t(H1) + W10 %*% mu11)
    W1[which(W1 <0)] <- W1

    T2 <- as.matrix((lambda2 / (lambda1 + eps)) * (A %*% W1))
    T2[which(T2 < 0)] <- 0
    numer <- t(W20) %*% E
    H2 <- H20 * (numer)/(t(W20) %*% W20 %*% H20 + eps)
    H2[which(H2 < 0 )] <- 0
    numer <- E%*%t(H2)
    mu21 <- diag(apply(W20*(numer + T2), MARGIN = 2, sum))
    mu22 <- diag(apply(W20*(W20 %*% (H2 %*% t(H2))), MARGIN = 2, sum))
    W2 <- W20 * (numer + T2 + W20 %*% mu22) / (eps + W20 %*% H2 %*% t(H2) + W20 %*% mu21)
    W2[which(W2 < 0)] <- 0

    dnorm1 <- norm((O-W1 %*% H1), type <- 'F')**2 + lambda1 * (norm(E - W2 %*% H2)**2)
    dw1 <- max(abs(W1-W10))/(max(abs(W10)) + eps)
    dh1 <- max(abs(H1-H10))/(max(abs(H10)) + eps)
    dw2 <- max(abs(W2-W20))/(max(abs(W20)) + eps)
    dh2 <- max(abs(H2-H20))/(max(abs(H20)) + eps)
    delta <- max(c(dw1, dh1, dw2, dh2))
    if (it>1){
      if(delta <= tolx){
        print(paste0('delta=', delta, " is smaller"))
        break
      }
      else if (abs(dnorm1 - dnorm0) < 0.01){
        print(paste0('Stop at the optimal solution'))
        break
        }
      else if (it == maxiter){
        break
      }
    }

    W10 <- W1
    W20 <- W2
    H10 <- H1
    H20 <- H2


    if (it %% 50 == 0 & nassign == 0){
      # print('start')
      S10 <- apply(H10, 2, function(x){which(x == max(x))})
      S20 <- apply(H20, 2, function(x){which(x == max(x))})
      FC1 <- matrix(0, nrow <- ngenes, ncol <- ncluster)
      FC2 <- matrix(0, nrow <- ngenes, ncol <- ncluster)
      # print('before')
      for (i in 1:ncluster){
        FC1[, i] <- apply(O, 1, function(x){sum(x[which(S10 == i)] > 2)}) / apply(O, 1, function(x){sum(x[which(S10 != i)] > 2)} + 1) * (sum(S10 == i)/sum(S10 != i) + 1)
        FC2[, i] <- apply(E, 1, function(x){sum(x[which(S20 == i)] > 0)}) / apply(E, 1, function(x){sum(x[which(S20 != i)] >0)} + 1) * (sum(S20 == i)/sum(S20 != i) + 1)
      }
      WP1 <- preprocessCore::normalize.quantiles.robust(FC1, use.median = TRUE)
      WP2 <- preprocessCore::normalize.quantiles.robust(FC2, use.median = TRUE)
      # print('after')
      S <- t(WP2) %*% A %*% WP1
      # print("run S")
      S <- 1/(1 + ((S-min(S)))/(max(S) - min(S)))
      err <- (dnorm0-dnorm1)/dnorm0

      print(paste0("Interation: ", it))
      # print(paste0('dnorm0 ',  round(dnorm0, 2)))
      # print(paste0('dnorm1 ', round(dnorm1, 2)))
      # print(paste0("error ", (dnorm0-dnorm1)/dnorm0))
      #
      nassign <- nassign + 1

      eq1 <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2)
      eq2 <- lambda1/2 * (norm(E - W2 %*% H2)**2)
      eq3 <- lambda2 * sum(diag(t(W2) %*% A %*% W1))
      score1<- eq1 + eq2 - eq3

      tmp_it <- c(tmp_it, it)
      tmp_score <- c(tmp_score, score1)

      print(paste0('objective function ', round(score1, 2)))

    }
    dnorm0 <- dnorm1
  }
  W1 <- W10
  H1 <- H10
  W2 <- W20
  H2 <- H20

  # Score <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2) + lambda1/2 * (norm(E - W2 %*% H2)**2) - lambda2 * sum(diag(t(W2) %*% A %*% W1))
  S10 <- apply(H1, 2, function(x){which(x == max(x))})
  S20 <- apply(H2, 2, function(x){which(x == max(x))})

  #### Reassignment
  FC1 <- matrix(0, nrow <- ngenes, ncol <- ncluster)
  FC2 <- matrix(0, nrow <- ngenes, ncol <- ncluster)
  for (i in 1:ncluster){
    FC1[, i] <- apply(O, 1, function(x){mean(x[which(S10 == i)] )})
    FC2[, i] <- apply(E, 1, function(x){mean(x[which(S20 == i)] )})
  }

  ### Calculate the cosine similarity
  similarity <- matrix(0, nrow = ncluster, ncol = ncluster)
  for (i in 1:ncluster){
    for (j in 1:ncluster){
      similarity[i, j] <- round(lsa::cosine(FC1[, i], FC2[, j]), 3)
    }
  }

  S = 1 - similarity
  assignment <- HungarianSolver(S)$pairs
  W2 <- W2[, assignment[, 2]]
  H2 <- H2[assignment[, 2], ]
  S1 <- apply(H1, 2, function(x){which(x == max(x))})
  S2 <- apply(H2, 2, function(x){which(x == max(x))})

  eq1 <- 0.5*(norm((O - W1 %*% H1), type <- 'F')**2)
  eq2 <- lambda1/2 * (norm(E - W2 %*% H2)**2)
  eq3 <- lambda2 * sum(diag(t(W2) %*% A %*% W1))
  score <- eq1 + eq2 - eq3
  print(paste0('eq1 ', round(eq1)))
  print(paste0('eq2 ', round(eq2)))
  print(paste0('eq3 ', round(eq3)))
  out <- paste("Iteration", c(as.character(it), as.character(eq1), as.character(eq2), as.character(eq3)))
  writeLines(out, txtfile)
  close(txtfile)
  end <- proc.time()
  print(paste0('Run time:', (end-start)[3][[1]], 'seconds'))

  write.csv(data.frame(iteration = tmp_it,
                       score = tmp_score), file = 'covengence_table.txt')

  return(list(H1, H2, W1, W2, S1, S2))
}

