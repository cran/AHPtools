## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(AHPtools)
library(readxl)

## -----------------------------------------------------------------------------
file <- system.file("extdata", "example_transport.xlsx", package = "AHPtools")

## -----------------------------------------------------------------------------
AHPstruc <- read_excel(file, sheet = "ahp")
print(AHPstruc, n=Inf)

## -----------------------------------------------------------------------------
pcm <- read_excel(file, sheet = "pcm")
as.data.frame(pcm)

## -----------------------------------------------------------------------------
tree <- viewAHPtree(AHPstruc)

## -----------------------------------------------------------------------------
print(tree, "level", limit = NULL)

## -----------------------------------------------------------------------------
w <- AHPweights(file, "ahp", "pcm")
alternatives_list <- lapply(w$AHPresult, function(x) x$alternatives)
alternatives_list

## -----------------------------------------------------------------------------
weights_list <- lapply(w$AHPresult, function(x) x$weights)
weights_list

## -----------------------------------------------------------------------------
file <- system.file("extdata", "example_automobile.xlsx", package = "AHPtools")

## -----------------------------------------------------------------------------
w <- AHPweights(file, "AHP", "PCM")
weights_list <- lapply(w$AHPresult, function(x) x$weights)
weights_list

## -----------------------------------------------------------------------------
alternatives_list <- lapply(w$AHPresult, function(x) x$alternatives)
alternatives_list

