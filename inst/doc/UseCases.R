## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# if (!requireNamespace("AHPtools", quietly = TRUE)) {
#     install.packages("AHPtools")
# }
library(AHPtools)

#devtools::load_all(".")
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(knitr)))
suppressWarnings(suppressMessages(library(kableExtra)))
library(knitr)
library(kableExtra)
library(tidyverse)
library(dplyr)

## ----echo=FALSE---------------------------------------------------------------
# Load the package namespace
ns <- asNamespace("AHPtools")

# Get the exported functions
exported <- getNamespaceExports("AHPtools")

# Get the function signatures
exported_functions <- mget(exported, envir = ns, ifnotfound = NA)
exported_functions <- exported_functions[sapply(exported_functions, is.function)]

# Print the function names and signatures
fnnm <- c(); fnsg <- c()
for (func_name in names(exported_functions)) {
  fnnm <- c(fnnm, func_name)
  fnsg <- c(fnsg, paste(deparse(args(exported_functions[[func_name]])),collapse=""))
}
# for (func_name in names(exported_functions)) {
#   cat(func_name, ":\t")
#   cat(str(args(exported_functions[[func_name]])))
# }
tbl <- data.frame("Function_Name"=fnnm, "Function_Signature"=fnsg)
kable(tbl)

## -----------------------------------------------------------------------------
help("CR", package="AHPtools")

## ----message=F, warning=F, exp1-----------------------------------------------
runs <- rep(3:12, each=100)
R <- unlist(lapply(runs, function(x) CR(createRandomPCM(x))$CR))
L <- lapply(runs, function(x) CR(createLogicalPCM(x)))
Lcr <- unlist(lapply(L, function(x) x$CR))
Lcons <- unlist(lapply(L, function(x) x$CRconsistent))
exp1DF <- data.frame(Order=runs,Random.PCM=R,Logical.PCM=Lcr,Inconsistent=!Lcons)

## ----message=FALSE------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
summaryExp1 <- exp1DF %>% group_by(Order) %>% 
  summarise(Random.PCM=mean(Random.PCM), Logical.PCM=mean(Logical.PCM)*100,
            "Logical PCMs"=mean(Inconsistent))
summaryExp1$Random.PCM <- round(summaryExp1$Random.PCM,3)
summaryExp1$Logical.PCM <- round(summaryExp1$Logical.PCM,3)
summaryExp1$`Logical PCMs` <- round(summaryExp1$`Logical PCMs`,2) 

kable(summaryExp1) %>%
  add_header_above(c(" "=1, "Average CR of 100"=2, "% Inconsistent"=1), line=FALSE)

## ----helper function----------------------------------------------------------
count_inversions <- function(x, y) {
  rx <- rank(x)
  ry <- rank(y)
  n <- length(rx)
  inversions <- 0

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Check if pair order differs between the two vectors
      if ((rx[i] - rx[j]) * (ry[i] - ry[j]) < 0) {
        inversions <- inversions + 1
      }
    }
  }

  return(inversions)
}

## -----------------------------------------------------------------------------
set.seed(93)
ind <- type <- ConsistencyRatio <- inverts <- vecRanks <- c()
oCR <- iCR <- cinv <- c()
for (i in 1:10) {
  p1 <- createRandomPCM(9)
  imp <- improveCR(p1)
  ind <- c(ind, i)
  #type <- c(type, c("Original", "Improved", ""))
  oCR <- c(oCR, round(imp$CR.original,4))
  iCR <- c(iCR, round(imp$suggestedCR,4))  
  spcm <- abs(Re(eigen(imp$suggestedPCM)$vectors[,1]))
  opcm <- abs(Re(eigen(p1)$vectors[,1]))
  cinv <- c(cinv, count_inversions(opcm, spcm))
  # cat(c("Iteration", i, ": ", ConsistencyRatio, "Inversions=", count_inversions(opcm, spcm), "\n"))
  # inverts <- c(inverts, c(imp$inversions, " ", ""))  
  # vecRanks <- c(vecRanks, c(imp$oriRank, imp$impRank, ""))
}
df <- data.frame(ind, oCR, iCR, cinv)
kable(df, linesep=FALSE) %>% kable_styling(position = "center")  

## -----------------------------------------------------------------------------
runs <- rep(c(5,7,9), each=20)
R <- unlist(lapply(runs, function(x) sensitivity(createRandomPCM(x))))
L <- unlist(lapply(runs, function(x) sensitivity(createLogicalPCM(x))))
exp2DF <- data.frame(Order=runs,Random.PCM=R,Logical.PCM=L)
exp2DF %>% group_by(Order) %>% summarise(R=mean(Random.PCM), L=mean(Logical.PCM), Rs=sd(Random.PCM), Ls=sd(Logical.PCM))


## -----------------------------------------------------------------------------
  pcmVec <- c(1/3,1/2,1/8,  1,1/6, 3)
  pcm <- createPCM(pcmVec)
  colnames(pcm) <- rownames(pcm) <- c('a1', 'a2', 'a3', 'a4')
  pcm

## -----------------------------------------------------------------------------
  #trdf <- triadReversal(pcm)
  trdf <- consEval(pcm)$triadsData
  trdf

## -----------------------------------------------------------------------------
  pcm3Vec <- c(1,3, 2)
  pcm3 <- createPCM(pcm3Vec)
  colnames(pcm3) <- rownames(pcm3) <- c('a2', 'a3', 'a1')
  pcm3

## -----------------------------------------------------------------------------
  max(0.4301402/1.1447143, 1.1447143/0.4301402)

## -----------------------------------------------------------------------------
trdf[,4]

## -----------------------------------------------------------------------------
consEval(pcm)

