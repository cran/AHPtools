eps <- 0.0001
consistentPCM <- matrix(c(1,9,1,1,1,  1/9,1,1/9,1/7,1/6,   1,9,1,1,1,
                          1,7,1,1,1,   1,6,1,1,1), nrow=5, byrow=TRUE)
rin <- c(0,0,0.52,0.89,1.11,1.25,1.35,1.40,1.45,1.49,1.52,1.54,1.56,1.58,1.59)
rth <- c(0,0,0.05,0.09,rep(0.1,11))
fs <- c(1/(9:1),2:9)
lim <- 500
lc2 <- lim*(lim-1)/2

notPCM <- function(PCM) {
  if (!setequal(diag(PCM),rep(1,nrow(PCM)))) return(TRUE)
  for (i in 1:(nrow(PCM)-1))
    for (j in (i+1):ncol(PCM))
      if (PCM[i,j]!=1/PCM[j,i]) return(TRUE)
  return(FALSE)
}

bestM <- function(pcm) {
  tSc <- c(1/(9:1),2:9)
  p <- pcm
  o <- nrow(pcm)
  bestMatrix <- diag(o)
  ep <- abs(Re(eigen(p)$vectors[,1]))
  for (r in 1:(o-1))
    for (c in (r+1):o) {
      b <- tSc[which.min(abs(ep[r]/ep[c]-tSc))[1]]
      bestMatrix[r, c] <- b
      bestMatrix[c, r] <- 1/b
    }
  return(bestMatrix)
}

randomPert <- function(val) {
  r <- which(rank(abs(val-fs))<=5)
  randomChoice <- fs[sample(r,1)]
  return(randomChoice)
}

perturb <- function(PCM) {
  pertPCM <- diag(rep(1,nrow(PCM)))
  for (i in 1:(nrow(PCM)-1))
    for (j in (i+1):nrow(PCM)) {
      r <- randomPert(PCM[i,j])
      pertPCM[i,j] <- r
      pertPCM[j,i] <- 1/r
    }
  return(pertPCM)
}


#' @title Saaty CR Consistency
#'
#' @description Computes and returns the Consistency Ratio (CR) for an input
#' PCM and its boolean status of CR consistency
#' @param PCM A pairwise comparison matrix
#' @returns A list of 3 elements, a boolean for the CR consistency of the
#' input PCM, the CR consistency value and the principal eigenvector
#' @importFrom stats runif
#' @examples
#' CR.pcm1 <- CR(matrix(
#'                  c(1,1,7,1,1, 1,1,5,1,1/3, 1/7,1/5,1,1/7,1/8, 1,1,7,1,1,
#'                  1,3,8,1,1), nrow=5, byrow=TRUE))
#' CR.pcm1
#' CR.pcm2 <- CR(matrix(
#'                   c(1,1/4,1/4,7,1/5, 4,1,1,9,1/4, 4,1,1,8,1/4,
#'                   1/7,1/9,1/8,1,1/9, 5,4,4,9,1), nrow=5, byrow=TRUE))
#' CR.pcm2
#' @export
CR <- function(PCM) {
  if (!is.matrix(PCM)) stop("Input is not a matrix")
  if (nrow(PCM)!=ncol(PCM)) stop("Input is not a square matrix")
  if (notPCM(PCM)) stop("Input is not a positive reciprocal matrix")
  CR <- ((Re(eigen(PCM)$values[1])-nrow(PCM))/(nrow(PCM)-1))/rin[nrow(PCM)]
  CR <- ifelse(abs(CR)<eps,0,CR)
  CRcons <- ifelse(CR<rth[nrow(PCM)],TRUE,FALSE)
  ev <- Re(eigen(PCM)$vectors[,1])
  return(list(CRconsistent=CRcons, CR=CR, eVec=ev))
}

#' @title Improve a PCM's consistency
#'
#' @description For an input pairwise comparison matrix (PCM) that is
#' inconsistent, this function returns a consistent PCM if possible,
#' with the relative preference for its alternatives as close as
#' possible to the original preferences, as in the principal right eigenvector.
#' @param PCM A pairwise comparison matrix
#' @returns A list of 4 elements, suggested PCM, a boolean for the CR
#' consistency of the input PCM, the CR consistency value, a boolean for the
#' CR consistency of the suggested PCM, the CR consistency value of the
#' suggested PCM
#' @importFrom stats runif
#' @examples
#' CR.suggest2 <- improveCR(matrix(
#'                  c(1,1/4,1/4,7,1/5, 4,1,1,9,1/4, 4,1,1,8,1/4,
#'                  1/7,1/9,1/8,1,1/9, 5,4,4,9,1), nrow=5, byrow=TRUE))
#' CR.suggest2
#' CR.suggest3 <- improveCR(matrix(
#'                  c(1,7,1,9,8, 1/7,1,1/6,7,9, 1,6,1,9,9, 1/9,1/7,1/9,1,5,
#'                  1/8,1/9,1/9,1/5,1), nrow=5, byrow=TRUE))
#' CR.suggest3
#' @export
improveCR <- function(PCM) {
  if (!is.matrix(PCM)) stop("Input is not a matrix")
  if (nrow(PCM)!=ncol(PCM)) stop("Input is not a square matrix")
  if (notPCM(PCM)) stop("Input is not a positive reciprocal matrix")
  CR <- ((Re(eigen(PCM)$values[1])-nrow(PCM))/(nrow(PCM)-1))/rin[nrow(PCM)]
  CR <- ifelse(abs(CR)<eps,0,CR)
  CRcons <- ifelse(CR<rth[nrow(PCM)],TRUE,FALSE)
  if (CRcons) stop("Input PCM is already CR consistent")
  sPCM <- bestM(PCM)
  sCR <- ((Re(eigen(sPCM)$values[1])-nrow(sPCM))/(nrow(sPCM)-1))/rin[nrow(sPCM)]
  sCR <- ifelse(abs(sCR)<eps,0,sCR)
  if(sCR > rin[nrow(sPCM)])
    stop("Input PCM though not CR consistent cannot be improved")
  sCRcons <- ifelse(sCR<rth[nrow(sPCM)],TRUE,FALSE)

  return(list(suggestedPCM=sPCM, CR.originalConsistency=CRcons,
              CR.original=CR, suggestedCRconsistent=sCRcons, suggestedCR=sCR))
}

#' @title Compute Sensitivity
#'
#' @description This function returns a sensitivity measure for an input
#' pairwise comparison matrix (PCM). Sensitivity is measured by Monte Carlo
#' simulation of 500 PCMs which are perturbations of the input PCM. The
#' perturbation algorithm makes a random choice from one of the 5 closest
#' items in the Fundamental Scale \{1/9, 1/8, ..... 1/2, 1, 2, ..... 8, 9\}
#' for each element in the PCM, ensuring the the pairwise reciprocity is
#' maintained. The sensitivity measure is the average Spearman's rank
#' correlation of the vector of ranks of the principal eigenvectors of
#' (i) the input PCM and (ii) the perturbed PCM. The average of the 500 such
#' rank correlations is reported as the measure of sensitivity.
#' @param PCM A pairwise comparison matrix
#' @returns The average Spearman's rank correlation between the principal
#' eigenvectors of the input and the perturbed PCMs
#' @importFrom stats runif
#' @examples
#' sensitivity1 <- sensitivity(matrix(
#'                  c(1,1/4,1/4,7,1/5, 4,1,1,9,1/4, 4,1,1,8,1/4,
#'                  1/7,1/9,1/8,1,1/9, 5,4,4,9,1), nrow=5, byrow=TRUE))
#' sensitivity1
#' sensitivity2 <- sensitivity(matrix(
#'                  c(1,7,1,9,8, 1/7,1,1/6,7,9, 1,6,1,9,9, 1/9,1/7,1/9,1,5,
#'                  1/8,1/9,1/9,1/5,1), nrow=5, byrow=TRUE))
#' sensitivity2
#' @export
sensitivity <- function(PCM) {
  if (!is.matrix(PCM)) stop("Input is not a matrix")
  if (nrow(PCM)!=ncol(PCM)) stop("Input is not a square matrix")
  if (notPCM(PCM)) stop("Input is not a positive reciprocal matrix")
  ev0 <- abs(Re(eigen(PCM)$vectors[,1]))
  d0 <- rank(-ev0)
  cs <- 0
  for (i in 1:lim) {
    c <- perturb(PCM)
    ev <- abs(Re(eigen(c)$vectors[,1]))
    d <- rank(-ev)
    cs <- cs + stats::cor(d0, d, method="spearman")
  }
  meanCor <- cs / lim
  return(meanCor)
}
