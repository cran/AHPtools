eps <- 0.0001
rin <- c(0,0,0.52,0.89,1.11,1.25,1.35,1.40,1.45,1.49,1.52,1.54,1.56,1.58,1.59)
rth <- c(0,0,0.05,0.09,rep(0.1,11))
fs <- c(1/(9:1),2:9)
# added 14.8.2023
fs2 <- sort(unique(as.vector(outer(fs, fs, "/"))))
lim <- 500
lc2 <- lim*(lim-1)/2
ord <- rep(3:12,each=2)
Q1 <- c(0.25,1.442,0.243,1.414,0.219,1.319,0.201,1.319,0.186,1.313,0.17,1.3,
        0.158,1.29,0.146,1.284,0.135,1.283,0.127,1.281)
Q2 <- c(0.376,1.747,0.325,1.348,0.281,1.341,0.247,1.314,0.222,1.304,0.2,1.294,
        0.182,1.297,0.169,1.303,0.155,1.296,0.144,1.292)
Q3 <- c(0.499,1.26,0.402,1.357,0.342,1.345,0.297,1.338,0.261,1.33,0.233,1.321,
        0.211,1.324,0.194,1.311,0.176,1.306,0.163,1.307)
Q4 <- c(0.746,1.4,0.691,1.402,0.627,1.381,0.563,1.353,0.482,1.339,0.435,1.34,
        0.424,1.331,0.331,1.326,0.293,1.322,0.278,1.319)
type <- rep(c('0D','1C'),10)
ec <- data.frame(cbind(ord,Q1,Q2,Q3,Q4,type))
logitCoefficients <- c(5.073699, 4.320408, -113.395213, -5.228240)
logitModel <- matrix(logitCoefficients)

notPCM <- function(PCM) {
  if (!setequal(diag(PCM),rep(1,nrow(PCM)))) return(TRUE)
  for (i in 1:(nrow(PCM)-1))
    for (j in (i+1):ncol(PCM))
      if (PCM[i,j]!=1/PCM[j,i]) return(TRUE)
  return(FALSE)
}

bestM <- function(pcm, granularityLow=TRUE) {
  if (granularityLow==TRUE) {
    opt <- fs
  }  else {
    opt <- fs2
  }
  #tSc <- c(1/(9:1),2:9)
  p <- pcm
  o <- nrow(pcm)
  bestMatrix <- diag(o)
  ep <- abs(Re(eigen(p)$vectors[,1]))
  for (r in 1:(o-1))
    for (c in (r+1):o) {
      b <- opt[which.min(abs(ep[r]/ep[c]-opt))[1]]
      bestMatrix[r, c] <- b
      bestMatrix[c, r] <- 1/b
    }
  return(bestMatrix)
}

randomPert <- function(val, granularityLow) {
  if (granularityLow==TRUE) {
    opt <- fs
  }  else {
    opt <- fs2
  }
  r <- which(rank(abs(val-opt))<=5)
  randomChoice <- opt[sample(r,1)]
  return(randomChoice)
}

perturb <- function(PCM, granularityLow=TRUE) {
  pertPCM <- diag(rep(1,nrow(PCM)))
  for (i in 1:(nrow(PCM)-1))
    for (j in (i+1):nrow(PCM)) {
      r <- randomPert(PCM[i,j], granularityLow)
      pertPCM[i,j] <- r
      pertPCM[j,i] <- 1/r
    }
  return(pertPCM)
}

mDev <- function(pcm, ppcm) {
  o <- nrow(pcm)
  gm <- 1
  for (r in 1:(o-1))
    for (c in (r+1):o) {
      rat <- ppcm[r,c] / pcm[r,c]
      rat <- ifelse(rat < 1, 1/rat, rat)
      gm <- gm * rat
    }
  mgm <- gm^(2/(o*(o-1)))
  return(inconGM=mgm)
}

#' @title Create a Random Pairwise Comparison Matrix of given order n
#' @description Create a Random Pairwise Comparison Matrix of order n for Analytic Hierarchy
#' Process using Saaty's 17-point Fundamental Scale for comparison ratios of the 
#' upper triangular elements of the order n PCM
#'
#' @param PCMsize The order of the random 'PCM' to be generated
#'
#' @returns A Pairwise Comparison Matrix corresponding to the given order
#' @importFrom stats runif
#' @examples
#' PCM <- createPCM(5);
#' PCM <- createPCM(11);
#' @export
createRandomPCM <- function(PCMsize) {
    n <- PCMsize
    fscale <- c(1/(9:2),1:9)
    vec <- sample(fscale, choose(n,2), replace=T)
    
    if (n!=as.integer(n)) {
        return(1)
    } else if (n<=2) {
        return(2)
    } else {
        pcm <- diag(n)
        vecPtr <- 0
        for (r in 1:(n-1))
            for (c in (r+1):n) {
                vecPtr <- vecPtr+1
                pcm[r,c] <- vec[vecPtr]
                pcm[c,r] <- 1/vec[vecPtr]
            }
    }
    return(pcm)
}

#' @param vec
#'
#' @title Create a Pairwise Comparison Matrix of order n for Analytic Hierarchy
#' Process from a vector of length n(n-1)/2 comparison ratios
#' @description Create a Pairwise Comparison Matrix of order n from a vector of
#' length n(n-1)/2 independent upper triangular elements
#'
#' @param vec The preference vector of length as the order of the 'PCM'
#'
#' @returns A Pairwise Comparison Matrix corresponding to the upper triangular
#' elements
#' @importFrom stats runif
#' @examples
#' PCM <- createPCM(c(1,2,0.5,3,0.5,2));
#' PCM <- createPCM(c(1,.5,2,1/3,4,2,.25,1/3,.5,1,.2,6,2,3,1/3));
#' @export
createPCM <- function(vec) {
    n <- (1+sqrt(1+8*length(vec)))/2
    if (n!=as.integer(n)) {
        return(1)
    } else {
        pcm <- diag(n)
        vecPtr <- 0
        for (r in 1:(n-1))
            for (c in (r+1):n) {
                vecPtr <- vecPtr+1
                pcm[r,c] <- vec[vecPtr]
                pcm[c,r] <- 1/vec[vecPtr]
            }
    }
    return(pcm)
}

#' @title Simulated Logical Pairwise Comparison Matrix for the
#' Analytic Hierarchy Process
#'
#' @description Creates a logical pairwise comparison matrix for the Analytic
#' Hierarchy Process such as would be created by a rational decision maker
#' based on a relative vector of preferences for the alternatives involved.
#' Choices of the pairwise comparison ratios are from the Fundamental Scale
#' and simulate a reasonable degree of error. The algorithm is modified from
#' a paper by Bose, A [2022], \doi{https://doi.org/10.1002/mcda.1784}
#'
#' @param ord The desired order of the Pairwise Comparison Matrix
#' @param prefVec The preference vector of length as the order of the
#' input matrix
#' @param granularityLow The Scale for pairwise comparisons; default (TRUE)
#' is the fundamental scale; else uses a more find grained scale, derived
#' from pairwise ratios of the elements of the Fundamental Scale.
#' @returns A Logical Pairwise Comparison Matrix
#' @importFrom stats runif
#' @examples
#' lPCM <- createLogicalPCM(3,c(1,2,3));
#' lPCM <- createLogicalPCM(5,c(0.25,0.4,0.1,0.05,0.2));
#' @export
createLogicalPCM <- function(ord, prefVec=rep(NA,ord), granularityLow=TRUE) {
  if (is.na(ord)) stop("The first parameter is mandatory")
  if (!is.numeric(ord) || ord %% 1 != 0)
    stop("The first parameter has to be an integer")
  if (!all(is.na(prefVec)) && !is.numeric(prefVec))
    stop("The second parameter has to be a numeric vector")
  if (!all(is.na(prefVec)) && length(prefVec)!=ord)
    stop("The length of the second parameter has to be the same as the first")

  if (granularityLow==TRUE) {
    opt <- fs
  }  else {
    opt <- fs2
  }
  # opt <- ifelse(granularityLow,fs,fs2)

  if (is.na(prefVec[1]))
    prefVec <- runif(ord)

  mperfect <- outer(prefVec, prefVec, "/")

  m <- bestM(mperfect,  granularityLow)
  # now creating a logical PCM
  for (r in 1:(ord-1)) {
    for (c in (r+1):ord) {
      m1 <- which.min(abs(opt-m[r,c]))
      m2 <- which.min(abs(opt[-m1]-m[r,c]))
      m3 <- which.min(abs(opt[-c(m1,m2)]-m[r,c]))
      # random choice from the nearest 3
      allChoices <- choices <- c(m1, m2, m3)
      if (m[r,c] >= 1) {
        choices <- allChoices[opt[allChoices] >= 1]
      } else if (m[r,c] < 1) {
        choices <- allChoices[opt[allChoices] <= 1]
      }
      m[r,c] <- sample(opt[choices],1)
      m[c,r] <- 1/m[r,c]
    }
  }
  return(logicalPCM=m)
}


#' @title Saaty CR Consistency
#'
#' @description Computes and returns the Consistency Ratio for an input
#' PCM and its boolean status of consistency based on Consistency Ratio
#'
#' @param typePCM boolean flag indicating if the first argument is a PCM or a
#' vector of upper triangular elements
#' @param PCM A pairwise comparison matrix
#'
#' @returns A list of 3 elements, a boolean for the 'CR' consistency of the
#' input 'PCM', the 'CR' consistency value and the principal eigenvector
#' @importFrom stats runif
#' @examples
#' CR.pcm1 <- CR(matrix(
#'                  c(1,1,7,1,1, 1,1,5,1,1/3, 1/7,1/5,1,1/7,1/8, 1,1,7,1,1,
#'                  1,3,8,1,1), nrow=5, byrow=TRUE))
#' CR.pcm1
#' CR.pcm1a <- CR(c(1,7,1,1, 5,1,1/3, 1/7,1/8, 1), typePCM=FALSE)
#' CR.pcm1a
#' CR.pcm2 <- CR(matrix(
#'                   c(1,1/4,1/4,7,1/5, 4,1,1,9,1/4, 4,1,1,8,1/4,
#'                   1/7,1/9,1/8,1,1/9, 5,4,4,9,1), nrow=5, byrow=TRUE))
#' CR.pcm2
#' CR.pcm2a <- CR(c(1/4,1/4,7,1/5, 1,9,1/4, 8,1/4, 1/9),typePCM=FALSE)
#' CR.pcm2a
#' @export
CR <- function(PCM,typePCM=TRUE) {
  if (!typePCM) {
    if (!is.vector(PCM)) stop("Input is not a vector of pairwise ratios")
    if (length(PCM)<3 | length(PCM)>66)
      stop("Input vector is not of appropriate length for a
           PCM of order 3 to 12")
    PCM <- createPCM(PCM)
    if (!is.matrix(PCM)) stop("Input vector does not have required values for
                              all upper triangular elements")
  } else {
    if (!is.matrix(PCM)) stop("Input is not a matrix")
    if (nrow(PCM)!=ncol(PCM)) stop("Input is not a square matrix")
    if (nrow(PCM)==2 | nrow(PCM)>12) stop("Input matrix should be
                                          of order 3 upto 12")
    if (notPCM(PCM)) stop("Input is not a positive reciprocal matrix")
  }
  CR <- ((Re(eigen(PCM)$values[1])-nrow(PCM))/(nrow(PCM)-1))/rin[nrow(PCM)]
  CR <- ifelse(abs(CR)<eps,0,CR)
  CRcons <- ifelse(CR<rth[nrow(PCM)],TRUE,FALSE)
  ev <- Re(eigen(PCM)$vectors[,1])
  return(list(CRconsistent=CRcons, CR=CR, eVec=ev))
}

#' @title Improve the CR consistency of a PCM
#'
#' @description For an input pairwise comparison matrix, PCM that is
#' inconsistent, this function returns a consistent PCM if possible,
#' with the relative preference for its alternatives as close as
#' possible to the original preferences, as in the principal right eigenvector.
#' @param PCM A pairwise comparison matrix
#' @param typePCM boolean flag indicating if the first argument is a PCM or a
#' vector of upper triangular elements
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
#' CR.suggest2a <- improveCR(c(1/4,1/4,7,1/5, 1,9,1/4, 8,1/4, 1/9),
#' typePCM=FALSE)
#' CR.suggest2a
#' CR.suggest3 <- improveCR(matrix(
#'                  c(1,7,1,9,8, 1/7,1,1/6,7,9, 1,6,1,9,9, 1/9,1/7,1/9,1,5,
#'                  1/8,1/9,1/9,1/5,1), nrow=5, byrow=TRUE))
#' CR.suggest3
#' @export
improveCR <- function(PCM,typePCM=TRUE) {
  if (!typePCM) {
    if (!is.vector(PCM)) stop("Input is not a vector of pairwise ratios")
    if (length(PCM)<3 | length(PCM)>66) stop("Input vector is not of
                                             appropriate length for a PCM of
                                             order 3 to 12")
    PCM <- createPCM(PCM)
    if (!is.matrix(PCM)) stop("Input vector does not have required values for
                              all upper triangular elements")
  } else {
    if (!is.matrix(PCM)) stop("Input is not a matrix")
    if (nrow(PCM)!=ncol(PCM)) stop("Input is not a square matrix")
    if (nrow(PCM)==2 | nrow(PCM)>12) stop("Input matrix should be of order
                                          3 upto 12")
    if (notPCM(PCM)) stop("Input is not a positive reciprocal matrix")
  }
  CR <- ((Re(eigen(PCM)$values[1])-nrow(PCM))/(nrow(PCM)-1))/rin[nrow(PCM)]
  CR <- ifelse(abs(CR)<eps,0,CR)
  CRcons <- ifelse(CR<rth[nrow(PCM)],TRUE,FALSE)
  #if (CRcons) stop("Input PCM is already CR consistent")
  sPCM <- bestM(PCM)
  sCR <- ((Re(eigen(sPCM)$values[1])-nrow(sPCM))/(nrow(sPCM)-1))/rin[nrow(sPCM)]
  sCR <- ifelse(abs(sCR)<eps,0,sCR)
  #if(sCR > rin[nrow(sPCM)])
  #  stop("Input PCM though not CR consistent cannot be improved")
  sCRcons <- ifelse(sCR<rth[nrow(sPCM)],TRUE,FALSE)

  return(list(suggestedPCM=sPCM, CR.originalConsistency=CRcons,
              CR.original=CR, suggestedCRconsistent=sCRcons, suggestedCR=sCR))
}

#' @title Compute Sensitivity
#'
#' @description This function returns a sensitivity measure for an input
#' pairwise comparison matrix, PCM. Sensitivity is measured by Monte Carlo
#' simulation of 500 PCMs which are perturbations of the input PCM. The
#' perturbation algorithm makes a random choice from one of the 5 closest
#' items in the Fundamental Scale \{1/9, 1/8, ..... 1/2, 1, 2, ..... 8, 9\}
#' for each element in the PCM, ensuring the the pairwise reciprocity is
#' maintained. The sensitivity measure is the average Spearman's rank
#' correlation of the vector of ranks of the principal eigenvectors of
#' (i) the input PCM and (ii) the perturbed PCM. The average of the 500 such
#' rank correlations is reported as the measure of sensitivity.
#' @param PCM A pairwise comparison matrix
#' @param typePCM boolean flag indicating if the first argument is a PCM or a
#' vector of upper triangular elements
#' @param granularityLow The Scale for pairwise comparisons; default (TRUE)
#' is the fundamental scale; else uses a more find grained scale, derived
#' from pairwise ratios of the elements of the Fundamental Scale.
#' @returns The average Spearman's rank correlation between the principal
#' eigenvectors of the input and the perturbed 'PCMs'
#' @importFrom stats runif
#' @examples
#' sensitivity2 <- sensitivity(matrix(
#'                  c(1,7,1,9,8, 1/7,1,1/6,7,9, 1,6,1,9,9, 1/9,1/7,1/9,1,5,
#'                  1/8,1/9,1/9,1/5,1), nrow=5, byrow=TRUE))
#' sensitivity2
#' @export
sensitivity <- function(PCM,typePCM=TRUE,granularityLow=TRUE) {
  if (!typePCM) {
    if (!is.vector(PCM)) stop("Input is not a vector of pairwise ratios")
    if (length(PCM)<3 | length(PCM)>66) stop("Input vector is not of
                                             appropriate length for a PCM of
                                             order 3 to 12")
    PCM <- createPCM(PCM)
    if (!is.matrix(PCM)) stop("Input vector does not have required values for
                              all upper triangular elements")
  } else {
    if (!is.matrix(PCM)) stop("Input is not a matrix")
    if (nrow(PCM)!=ncol(PCM)) stop("Input is not a square matrix")
    if (nrow(PCM)==2 | nrow(PCM)>12) stop("Input matrix should be of order
                                          3 upto 12")
    if (notPCM(PCM)) stop("Input is not a positive reciprocal matrix")
  }
  ev0 <- abs(Re(eigen(PCM)$vectors[,1]))
  d0 <- rank(-ev0)
  cs <- 0
  for (i in 1:lim) {
    c <- perturb(PCM, granularityLow)
    ev <- abs(Re(eigen(c)$vectors[,1]))
    d <- rank(-ev)
    cs <- cs + stats::cor(d0, d, method="spearman")
  }
  meanCor <- cs / lim
  return(meanCor)
}

#' @title Evaluate Revised Consistency
#'
#' @description This function returns the revised consistency classification
#' for a PCM, evaluated by comparison with the threshold of consistency for
#' intentional PCMs in the same preference heterogeneity quartile. The measure
#' for inconsistency is the geometric mean of ratios in comparison with the
#' corresponding benchmark PCM.
#'
#' @param PCM A pairwise comparison matrix
#' @param typePCM boolean flag indicating if the first argument is a PCM or a
#' vector of upper triangular elements
#' @returns A list of four elements,
#' revCons = the revised consistency classification,
#' inconGM = the Geometric Mean measure of inconsistency with the best 'PCM',
#' dQrtl = the preference heterogeneity quartile for the normalized
#' eigenvector, and diff = the preference heterogeneity measure
#' @importFrom stats runif
#' @examples
#' revCon1 <- revisedConsistency(matrix(
#'                  c(1,1/4,1/4,7,1/5, 4,1,1,9,1/4, 4,1,1,8,1/4,
#'                  1/7,1/9,1/8,1,1/9, 5,4,4,9,1), nrow=5, byrow=TRUE))
#' revCon1
#' revCon2 <- revisedConsistency(c(7,1,9,8, 1/6,7,9, 9,9, 5), typePCM=FALSE)
#' revCon2
#' @export
revisedConsistency <- function(PCM,typePCM=TRUE) {
  if (!typePCM) {
    if (!is.vector(PCM)) stop("Input is not a vector of pairwise ratios")
    if (length(PCM)<3 | length(PCM)>66)
      stop("Input vector is not of appropriate length for a PCM of
           order 3 to 12")
    PCM <- createPCM(PCM)
    if (!is.matrix(PCM))
      stop("Input vector does not have required values for all
           upper triangular elements")
  } else {
    if (!is.matrix(PCM)) stop("Input is not a matrix")
    if (nrow(PCM)!=ncol(PCM)) stop("Input is not a square matrix")
    if (nrow(PCM)==2 | nrow(PCM)>12)
      stop("Input matrix should be of order 3 upto 12")
    if (notPCM(PCM)) stop("Input is not a positive reciprocal matrix")
  }
  evector <- abs(Re(eigen(PCM)$vectors[,1]))
  evector <- evector/sum(evector)
  diff <- max(evector)[1] - min(evector)[1]
  d <- as.numeric(unname(ec[ec$ord==nrow(PCM) & ec$type=='0D',2:5]))
  inc <- ec[ec$ord==nrow(PCM) & ec$type=='1C',2:5]
  # The true part of this is statement added to solve
  # the problem of min(which(d>diff)) not having any argument
  if (max(d)[1]<diff) {
    dQrtl <- "Q4"
    inconGM <- mDev(PCM, bestM(PCM))
    inconsThreshold <- as.numeric(inc[1])
    revCons <- inconGM <= inconsThreshold
  } else {
    column <- min(which(d>diff))
    # Need to add the below line to take care of PCMs with eigenvalue = 0
    column <- ifelse(is.infinite(column),1,column)
    inconsThreshold <- as.numeric(inc[column])
    inconGM <- mDev(PCM, bestM(PCM))
    dQrtl <- paste0("Q",column)
    revCons <- inconGM <= inconsThreshold
  }
  return(list(revCons=revCons,inconGM=inconGM,dQrtl=dQrtl,diff=diff))
}

reversals <- function(arr1, arr2) {
  # arr1 : slice from parent
  # arr2 : child
  n <- length(arr1)
  l <- list()

  signRev <- list()
  vrev <- list()
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      r1 <- arr1[i]/arr1[j]; r2 <- arr2[i]/arr2[j]
      revTrue <- ((r1 > 1) & (r2 < 1)) |  ((r1 < 1) & (r2 > 1))  # Reversal 
      if (revTrue) {
         k <- setdiff(names(arr1), c(names(arr1[i]), names(arr1[j])))
         l <- append(l,list(c(arr1[i], arr1[j], arr1[k], arr2[i], arr2[j], arr2[k])))
         vrev <- append(vrev,max(r1/r2, r2/r1))
      }
      
    } # for j
  }   # for i
  return(list(vrev=vrev,rev=l))
}

#' @importFrom utils combn
triadReversal <- function(PCM) {
  df <- data.frame()
  rownames(PCM) <- colnames(PCM) <- paste0("a",1:nrow(PCM))
  ev <- abs(Re(eigen(PCM)$vectors[,1]))
  evl <- Re(eigen(PCM)$values[1])
  names(ev) <- colnames(PCM)
  ch <- combn(1:nrow(PCM), 3)
  for (i in 1:ncol(ch)) {
    submatrix <- PCM[ch[,i], ch[,i]]
    e <- Re(eigen(submatrix)$vectors[,1])
    names(e) <- colnames(submatrix)
    ev2 <- ev[names(e)]
    inv <- reversals(ev2, e)

    if (length(inv$vrev)>0) {
      for (k in 1:length(inv$vrev)) {
        df <- rbind(df, c(unlist(names(inv$rev[[k]]))[1:3],unlist(abs(inv$vrev[[k]])),unlist(abs(inv$rev[[k]]))))
      }   # for (k in 1:length(inv$vrev))
    }     # if (length(inv$vrev)>0)
  }       # for (i in 1:ncol(ch))
  if (nrow(df)>0) {
    colnames(df) <- c("triadE1", "triadE2", "triadE3", "prefRev", "pcmWeightE1", "pcmWeightE2", "pcmWeightE3", "triadWeightE1", "triadWeightE2", "triadWeightE3")
    for (numCol in 4:10) 
      df[,numCol] <- as.numeric(df[,numCol])
  } else
    df <- NULL
  return(df)
}

#' @title Find consistency of a PCM based on Preference Reversals
#'
#' @description This function finds all triad based preference reversals for a PCM.
#' #' Triads are subsets of 3 elements chosen from the 'n' alternatives of an order-n
#' PCM. A triad reversal is said to occur if any two elements of the order-3 PCM
#' show a reversal in preference with the corresponding elements of the full eigenvector.
#' 
#' This returns a list of values related to triad Preference Reversal consistency.
#' The fourth item of the list is a data frame of triads where reversals are seen,
#' with the logit Consistency probability measure based on them, the proportion of
#' reversals, the maximum reversal and the data frame with the details.
#' @param pcm A pairwise comparison matrix
#' @returns A list of four elements,
#' logitConsistency = the probability that the PCM is consistent,
#' prop3Rev = the proportion of triad-based preference reversals for the PCM,
#' max3Rev = the maximum triad-based preference reversal for the PCM,
#' triadsData = a data frame with 8 columns, providing the full data of preference reversals
#'        (1) triadE1 alternative 1  in the triad; e.g. a4 for the fourth alternative
#'        (2) triadE2 alternative 2 in the triad
#'        (3) triadE3 alternative 3 in the triad
#'        (4) pref3Rev measure of the intensity of preference reversal for the particular triad
#'        (5) pcmWeightE1 eigen weight of alternative triadE1 from the entire eigenvector
#'        (6) pcmWeightE2 eigen weight of alternative triadE2 from the entire eigenvector
#'        (7) triadWeightE1 eigen weight of alternative triadE1 from the order-3 sub matrix
#'        (8) triadWeightE2 eigen weight of alternative triadE2 from the order-3 sub matrix
#' @examples
#' pcm1 <- matrix(c(1,1,2,1,2,2, 1,1,1,1/3,1,1, 1/2,1,1,1,1,1, 1,3,1,1,2,1,
#'                  1/2,1,1,1/2,1,1/4, 1/2,1,1,1,4,1), nrow=6, byrow=TRUE)
#' cons1 <- consEval(pcm1)
#' cons1
#' pcm2 <- matrix(c(1,1/6,1/5,1/2,1/6,1/3,1/3,1/8,  6,1,1,3,1,1,2,1,
#'                  5,1,1,3,1/3,1,2,1/2, 2,1/3,1/3,1,1/5,1/2,1/2,1/4,  
#'                  6,1,3,5,1,1,2,1,  3,1,1,2,1,1,1,1/3,
#'                  3,1/2,1/2,2,1/2,1,1,1/4,  8,1,2,4,1,3,4,1), 
#'                  nrow=8, byrow=TRUE)
#' cons2 <- consEval(pcm2)
#' cons2
#' pcm3 <- createLogicalPCM(7) 
#' cons3 <- consEval(pcm3)
#' print(paste(formatC(cons3$logitConsistency,format="e",digits=4), 
#'          formatC(cons3$prop3Rev,format="f",digits=4), 
#'          formatC(cons3$max3Rev,format="f", digits=4)))
#' @export
consEval <- function(pcm) {
  a <- triadReversal(pcm)
  if (!is.null(a)) {
    b <- data.frame(1, order=nrow(pcm), prop3Rev=nrow(a)/(choose(nrow(pcm),3)*3), max3Rev=max(a$prefRev))
    # c <- predict(logitModel, newdata=b)
    c <- as.numeric(as.matrix(b,nrow=1) %*% matrix(logitModel))
    d <- 1/(1/exp(c) + 1)
    logitConsistent <- ifelse(d>0.5,TRUE, FALSE)
    triadsData <- a[,-c(7,10)]
  } else {
    d <- 1
    b <- c(1,order=nrow(pcm),0,1)
    triadsData <- NULL
  }
  #return(list(unname(d), unname(logitConsistent), unname(b[c2,3)]), a[,-c(7,10)]))
  return(list(logitConsistency=unname(d), prop3Rev=as.numeric(b[3]), max3Rev=as.numeric(b[4]), triadsData=triadsData))
}

#' @importFrom readxl read_excel
readAHP <- function(excelFile, ahpSheet=1, pcmSheet=2) {
  ahp <- read_excel(excelFile, sheet = ahpSheet)
  pcm <- read_excel(excelFile, sheet = pcmSheet) 
  return(list(ahp, pcm))
}

#' @keywords internal
#' @noRd
pcmCreate <- function(node, resp, ahpSh, pcmSh) {
  grepStr <- paste("^",node,"\\.\\d+$", sep="")
  utCols <- grep(grepStr, names(pcmSh))
  pcm <- createPCM(unlist(pcmSh[resp,utCols]))
  ev <- abs(Re(eigen(pcm)$vector[,1]))
  ev <- ev/sum(ev)
  names(ev) <- strsplit(unlist(ahpSh[ahpSh[,1]==node,3]), ",")[[1]]
  return(ev)
}

#' @keywords internal
#' @noRd
ahpStructure <- function(df.ahp) {
  ahp <- list()
  for (i in 1:nrow(df.ahp)) {
    chString = unname(unlist(df.ahp[i,"Children_Ordered"]))
    ahp[[i]] <- list(
      node = unname(unlist(df.ahp[i,"Node"])),
      parent = unname(unlist(df.ahp[i,"Parent"])),
      children = strsplit(chString,",")[[1]]
    )
  }
  return(ahp)
}

#' @keywords internal
#' @noRd
# Recursive function to build the tree
add_children <- function(current_node, current_name, ahp_list) {
  # Find the matching node in ahp
  node_entry <- Filter(function(x) x$node == current_name, ahp_list)
  
  # Safety: check if found
  if (length(node_entry) == 0) return()
  
  children <- node_entry[[1]]$children
  
  # If children is a comma-separated string, split it
  if (is.character(children) && length(children) == 1 && grepl(",", children)) {
    children <- strsplit(children, ",")[[1]]
  }
  
  # Add each child and recurse
  for (child_name in children) {
    child_node <- current_node$AddChild(child_name)
    add_children(child_node, child_name, ahp_list)
  }
}


#' Create AHP Tree Structure
#'
#' Builds a hierarchical tree from a flat AHP representation.
#'
#' @param ahp A data frame with the AHP structure including Node and Parent columns
#'
#' @returns A `Node` object (from the `data.tree` package) representing the full AHP tree
#'
#' @examples
#' file <- system.file("extdata", "example_transport.xlsx", package = "AHPtools")
#' AHPstruc <- readxl::read_excel(file, sheet = "ahp")
#' tree <- viewAHPtree(AHPstruc)
#' print(tree, "level", limit=NULL)
#' @details
#' For an overview and examples, please see the associated vignette:
#' `vignette("viewAHPtree", package = "AHPtools")`
#' @importFrom data.tree Node
#' @export
viewAHPtree <- function(ahp) {
  ahpTree <- ahpStructure(ahp)
  rootNode <- unname(unlist(ahp[is.na(ahp$Parent), "Node"]))
  root <- Node$new(rootNode)
  add_children(root, rootNode, ahpTree)
  return(root)
}


#' @keywords internal
#' @noRd
validateAHP <- function(ahp, pcm) {
  cols <- c()
  for (i in 1:nrow(ahp)) {
    node <- unlist(ahp[i,1])
    children <- strsplit(unlist(ahp[i,3]), ",")
    ut <- choose(length(unlist(children)),2)
    cols <- c(cols, paste(node,1:ut, sep="."))
  }
  return(cols)
}

#' @title Compute weights for Alternatives and lowest level sub criteria in AHP responses
#'
#' @description This function reads an Excel file with two required Sheets, viz. 
#' Sheet 1:  for the AHP structure, with three columns as follows:
#' Sheet 2:  for the upper triangular elements of the PCMs that are part of the AHP hierarchy                          
#' 
#' This returns a list of two values:
#'  (1) a printable AHP tree excluding the alternatives, if any
#'  (2) the list of weights for the lowest level subcriteria, and weights of alternatives 
#'            if exists
#' @param ExcelPath for the Excel file containing the AHP structure and the required PCMs
#' @param AHPsheet  for the AHP structure, with three required columns, viz.
#' Column 1: Node:             the node names for all nodes that have child nodes 
#' Column 2: Parent:           the parent node for the Node in Column 1
#' Column 3: Children_Ordered: the child nodes for the Node in Column 1.
#'                             these are comma separated strings, and correspond to the
#'                             ordered upper triangular elements of the PCM in Sheet 2
#' @param PCMsheet  for the PCMs that are part of the AHP. The upper triangular matrix elements
#'                  are provided for each PCM, so that a nxn PCM has n(n-1)/2 entries.
#'                  These entries have column names starting with the AHP node name with
#'                  respect to which the child elements are being compared, followed by a 
#'                  dot (.) and a sequence of numbers from 1 to n(n-1)/2 for the PCM elements
#' @returns A list of two items,
#' (i)  AHPtree     which is a printable tree object constructed from the user-provided AHP structure
#' (ii) AHPresult   the list of weights for the lowest level subcriteria, and weights of alternatives 
#'                  if exists
#' @examples
#' file <- system.file("extdata", "example_transport.xlsx", package = "AHPtools")
#' results <- AHPweights(file, "ahp", "pcm")
#' print(results)
#' @details
#' For an overview and examples, please see the associated vignette:
#' `vignette("AHPweights", package = "AHPtools")`
#' @export
AHPweights <- function(ExcelPath, AHPsheet, PCMsheet) {   

  ra <- readAHP(excelFile= ExcelPath, 
                ahpSheet=AHPsheet, pcmSheet =PCMsheet)

  cols <- validateAHP(ra[[1]], ra[[2]])
  if (!identical(sort(cols),sort(colnames(ra[[2]])))) {
    stop(paste("Column names mis-specified in header of ", PCMsheet, sep=""))
    #valid <- FALSE
  } 

    ct <- viewAHPtree(ra[[1]])
    
    AHP_result <- list()
    for (rsp in 1:nrow(ra[[2]]))  {
      ra[[1]]$evec <- NA
      ra[[1]]$lWeight <- NA
      ra[[1]][is.na(ra[[1]]$Parent),"lWeight"] <- 1
        for (node in unlist(ra[[1]]$Node)) {
          pg <- pcmCreate(node,rsp, ra[[1]], ra[[2]])
          parent <- unlist(ra[[1]][ra[[1]]$Node == names(pg)[1], "Parent"]) 
          if (length(parent)>0) {
            ra[[1]][ra[[1]]$Node == parent,]$evec <- list(pg)
          } 
          if (! names(pg)[1] %in% unlist(ra[[1]]$Node)) {
            ra[[1]][ra[[1]]$Node == node,  ]$evec <- list(pg)
          }
          for (tnode in names(pg)) {
            ra[[1]][ra[[1]]$Node==tnode, "lWeight"] <- pg[tnode]
          }
        }
      
      ra[[1]]$GWeight <- ra[[1]]$lWeight 
      raFull <- ra[[1]]
      
      # Create Global Weights iteratively
      nparent <- c(unname(unlist(raFull[is.na(raFull$Parent), "Node"])))
      while (length(nparent)>0) {
        n2parent <- c()
        ra1 <- raFull[raFull$Parent %in% nparent,]
        for (i in 1:nrow(ra1)) {
          parent <- unlist(ra1[i,"Parent"])
          if (parent %in% nparent) {
              pweight <- raFull[raFull$Node==parent, "GWeight"]
              nweight <- ra1[i, "GWeight"]
              thisNode <- unlist(ra1[i,"Node"])
              raFull[raFull$Node==thisNode, "GWeight"] <- nweight * pweight
              n2parent <- c(n2parent, thisNode)
          }
        }
        nparent <- unique(ra[[1]]$Parent[ra[[1]]$Parent %in% n2parent])
      }
      
      level <- 0
      raFull[is.na(raFull$Parent),"level"] <- level
      lParents <- unique(unlist(raFull[is.na(raFull$Parent), "Node"]))
      while (length(lParents)>0) {
        level <- level + 1
        raFull[raFull$Parent %in% lParents, "level"] <- level  
        lParents <- unique(unlist(raFull[raFull$Parent %in% lParents,"Node"]))
      }
    
      raFull$weights <- NA
      for (i in 1:nrow(raFull)) {
        raFull$weights[i] <- list(raFull$GWeight[i] * unlist(raFull$evec[i]))
      }
    
      #flat <- unlist(raFull$weights[7:15])  # for Mathi data
      flat <- unlist(raFull$weights)
    
      # Sum by names
      result <- tapply(flat, names(flat), sum)
    
      # Which Nodes do not appear in the Parent list
      # These are the bottom level nodes in the AHP
      ch <- which(! raFull$Node %in% raFull$Parent)
      leaves <- raFull$Node[ch]
    
      # These are the alternatives
      alts <- c()
      for (i in 1:nrow(raFull)) {
        e <- names(unlist(raFull[i, "evec"]))
        cleaned <- sub("^evec\\.", "", e)
        alts <- c(alts, cleaned[! (cleaned %in% unlist(raFull$Node)) & ! (cleaned %in% alts)] )
      }  
    
      n <- unlist(raFull$Node)
      p <- unlist(raFull$Parent)
      r <- names(result)
      
      alts <- names(table(names(flat))[table(names(flat))>1]) 
      leaves <- setdiff(r,alts)
      leaves <- if (length(alts)==0) setdiff(leaves,n) else setdiff(leaves,p)
    
      AHP_result[[rsp]] <- 
         list(
            weights = result[leaves],
            alternatives = if (length(alts) > 0) result[alts] else NA
          )
    }  
    return(list(AHPtree=ct, AHPresult=AHP_result))
}