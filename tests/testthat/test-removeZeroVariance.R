library(geneExpr)
context("removeZeroVariance")

# create a test matrix
set.seed(378)
testMat <- matrix( sample(1:100,50), ncol=10)
# add some rows/columns with zero variance
testMat[ 3, ] <- rep(2,10)
testMat[ , 6 ] <- rep(2,5)
testMat[ , 8 ] <- rep(2,5)

test_that("check names", {
  expect_error( RemoveZeroVariance(testMat), 
                "The supplied matrix does not have row names!")
  expect_error( RemoveZeroVariance(testMat, cols=TRUE), 
                "The supplied matrix does not have column names!")
})
