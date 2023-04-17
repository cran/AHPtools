test_that("errors", {
  testthat::expect_error(
    CR(as.matrix(rep(1,9),byrow=TRUE,nrow=3)),
    "Input is not a square matrix"
  )

  testthat::expect_error(
    CR(matrix(rep(1.5,9),byrow=TRUE,nrow=3)),
    "Input is not a positive reciprocal matrix"
  )

  testthat::expect_equal(
    length(CR(matrix(rep(1,9),byrow=TRUE,nrow=3)))==3,
    TRUE
  )

  testthat::expect_equal(
    class(CR(matrix(rep(1,9),byrow=TRUE,nrow=3))[[1]])=="logical",
    TRUE
  )

  testthat::expect_equal(
    length(CR(matrix(rep(1,9),byrow=TRUE,nrow=3))[[3]])==3,
    TRUE
  )

  testthat::expect_equal(
    CR(matrix(rep(1,9),byrow=TRUE,nrow=3))[[2]]==0 &
      CR(matrix(rep(1,9),byrow=TRUE,nrow=3))[[1]],
    TRUE
  )

  testthat::expect_error(
    CR(c(1:16)),
    "Input is not a matrix"
  )

  testthat::expect_error(
    CR(matrix(c(1:24),nrow=4,byrow=TRUE)),
    "Input is not a square matrix"
  )

  testthat::expect_error(
    CR(matrix(c(1:36),nrow=6,byrow=TRUE)),
    "Input is not a positive reciprocal matrix"
  )

  testthat::expect_error(
    improveCR(c(1:16)),
    "Input is not a matrix"
  )

  testthat::expect_error(
    improveCR(matrix(c(1:24),nrow=4,byrow=TRUE)),
    "Input is not a square matrix"
  )

  testthat::expect_error(
    improveCR(matrix(c(1:36),nrow=6,byrow=TRUE)),
    "Input is not a positive reciprocal matrix"
  )

  testthat::expect_error(
    improveCR(matrix(c(1,9,1,1,1,  1/9,1,1/9,1/7,1/6,   1,9,1,1,1,
                       1,7,1,1,1,   1,6,1,1,1), nrow=5, byrow=TRUE)),
    "Input PCM is already CR consistent"
  )

  testthat::expect_equal(
    improveCR(matrix(c(1,9,1,1/4,1,  1/9,1,1/9,7,1/6,   1,9,1,1,1,
                       4,1/7,1,1,1,   1,6,1,1,1), nrow=5, byrow=TRUE))[[2]],
    FALSE
  )

  testthat::expect_equal(
    improveCR(matrix(c(1,9,1,1/4,1,  1/9,1,1/9,7,1/6,   1,9,1,1,1,
                       4,1/7,1,1,1,   1,6,1,1,1), nrow=5, byrow=TRUE))[[4]],
    TRUE
  )

  testthat::expect_equal(
    nrow(improveCR(matrix(c(1,9,1,1/4,1,  1/9,1,1/9,7,1/6,   1,9,1,1,1,
                            4,1/7,1,1,1,   1,6,1,1,1),
                          nrow=5, byrow=TRUE))[[1]]),
    5
  )

  testthat::expect_error(
    sensitivity(c(1:16)),
    "Input is not a matrix"
  )

  testthat::expect_error(
    sensitivity(matrix(c(1:24),nrow=4,byrow=TRUE)),
    "Input is not a square matrix"
  )

  testthat::expect_error(
    sensitivity(matrix(c(1:36),nrow=6,byrow=TRUE)),
    "Input is not a positive reciprocal matrix"
  )

  testthat::expect_equal(
    sensitivity(matrix(c(1,9,1,1/4,1,  1/9,1,1/9,7,1/6,   1,9,1,1,1,
                         4,1/7,1,1,1,   1,6,1,1,1), nrow=5, byrow=TRUE))>0,
    TRUE
  )

  testthat::expect_equal(
    sensitivity(matrix(c(1,9,1,1/4,1,  1/9,1,1/9,7,1/6,   1,9,1,1,1,
                         4,1/7,1,1,1,   1,6,1,1,1), nrow=5, byrow=TRUE))<1,
    TRUE
  )

})
