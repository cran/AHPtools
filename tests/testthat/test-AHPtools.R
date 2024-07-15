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

  testthat::expect_lte(
    sensitivity(matrix(c(1,9,1,1/4,1,  1/9,1,1/9,7,1/6,   1,9,1,1,1,
                         4,1/7,1,1,1,   1,6,1,1,1), nrow=5, byrow=TRUE)),
    1
  )

  testthat::expect_true(
    is.matrix(createPCM(c(9,1,1/4,1,  1/9,7,1/6,  1,1,  1)))
  )

  testthat::expect_equal(
    nrow(createPCM(1:6)),
    4
  )
  
  testthat::expect_error(
    createLogicalPCM(NA, c(1,2,3)),
    "The first parameter is mandatory"
  )

  testthat::expect_error(
    createLogicalPCM(3.25),
    "The first parameter has to be an integer"
  )

  testthat::expect_error(
    createLogicalPCM(3, c(1, "a", 2.5)),
    "The second parameter has to be a numeric vector"
  )

  testthat::expect_error(
    createLogicalPCM(3, c(1,2,.75, 3)),
    "The length of the second parameter has to be the same as the first"
  )
  
  testthat::expect_in(
    revisedConsistency(matrix(
                        c(1,1/4,1/4,7,1/5, 4,1,1,9,1/4, 4,1,1,8,1/4,
                        1/7,1/9,1/8,1,1/9, 5,4,4,9,1), nrow=5, byrow=TRUE))[1],
    c(TRUE, FALSE)
    )
  
  testthat::expect_in(
    revisedConsistency(c(1/4,1/4,7,1/5, 1,9,1/4, 8,1/4, 1/9),typePCM=FALSE)[1],
    c(TRUE, FALSE)
    )

  testthat::expect_in(
    revisedConsistency(c(9,1/9,3,1,1,1/4, 1/7,1/2,1/5,1/4,1/9, 6,5,2,1/9, 
                        1/7,1/6,1/9, 1/7,1/8, 1/8),typePCM=FALSE)[1],
    c(TRUE, FALSE)
    )
  testthat::expect_gt(
    consEval(matrix(c(1,1/6,1/5,1/2,1/6,1/3,1/3,1/8,  6,1,1,3,1,1,2,1,  5,1,1,3,1/3,1,2,1/2,
                 2,1/3,1/3,1,1/5,1/2,1/2,1/4,  6,1,3,5,1,1,2,1,  3,1,1,2,1,1,1,1/3,
                 3,1/2,1/2,2,1/2,1,1,1/4,  8,1,2,4,1,3,4,1), nrow=8, byrow=T))$logitConsistency,
    0.9999
    )
  testthat::expect_gte(
    consEval(createLogicalPCM(8))$logitConsistency,
    0
    )
  testthat::expect_lte(
    consEval(createLogicalPCM(5))$logitConsistency,
    1.00
    )
  testthat::expect_gte(
    consEval(createLogicalPCM(9))$max3Rev,
    1
    )
  testthat::expect_null(
    consEval(matrix(c(1,3,1/2,1,  1/3,1,1/7,1/4,  2,7,1,1,  1,4,1,1), 
        nrow=4, byrow=T))$triadsData
    )
  testthat::expect_equal(
    consEval(matrix(c(1,3,1/2,1,  1/3,1,1/7,1/4,  2,7,1,1,  1,4,1,1), 
      nrow=4, byrow=T))$prop3Rev,
    0
    )
  
})
