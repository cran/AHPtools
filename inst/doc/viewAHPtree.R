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
tree <- viewAHPtree(AHPstruc)

## -----------------------------------------------------------------------------
print(tree, "level", limit = NULL)

## -----------------------------------------------------------------------------
file <- system.file("extdata", "example_automobile.xlsx", package = "AHPtools")
AHPstruc2 <- read_excel(file, sheet = "AHP")
print(AHPstruc2, n=Inf)
tree <- viewAHPtree(AHPstruc2)
print(tree, "level", limit = NULL)

