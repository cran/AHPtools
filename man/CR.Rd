% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AHPtools.R
\name{CR}
\alias{CR}
\title{Saaty CR Consistency}
\usage{
CR(PCM)
}
\arguments{
\item{PCM}{A pairwise comparison matrix}
}
\value{
A list of 3 elements, a boolean for the CR consistency of the
input PCM, the CR consistency value and the principal eigenvector
}
\description{
Computes and returns the Consistency Ratio (CR) for an input
PCM and its boolean status of CR consistency
}
\examples{
CR.pcm1 <- CR(matrix(
                 c(1,1,7,1,1, 1,1,5,1,1/3, 1/7,1/5,1,1/7,1/8, 1,1,7,1,1,
                 1,3,8,1,1), nrow=5, byrow=TRUE))
CR.pcm1
CR.pcm2 <- CR(matrix(
                  c(1,1/4,1/4,7,1/5, 4,1,1,9,1/4, 4,1,1,8,1/4,
                  1/7,1/9,1/8,1,1/9, 5,4,4,9,1), nrow=5, byrow=TRUE))
CR.pcm2
}
