library(geneExpr)
context("removeZeroVariance")

# create a test matrix
set.seed(378)
testMat <- matrix( sample(1:100,50), ncol=10)
# add some rows/columns with zero variance
testMat[ 3, ] <- rep(2,10)
testMat[ , 6 ] <- rep(2,5)
testMat[ , 8 ] <- rep(2,5)

test_that("check matrix dim names", {
  expect_error( RemoveZeroVariance(testMat), 
                "The supplied matrix does not have row names!")
  rownames(testMat) <- paste0('geneId_', 1:5)
  expect_error( RemoveZeroVariance(testMat, cols=TRUE), 
                "The supplied matrix does not have column names!")
})

rownames(testMat) <- paste0('geneId_', 1:5)
colnames(testMat) <- paste0('sampleId_', 1:10)

test_that("check returned matrix: Rows only", {
  resultMat <- testMat[ c(1:2,4:5), ]
  # remove zero var, rows only
  res <- RemoveZeroVariance(testMat)
  expect_equal( nrow(res$matrix), nrow(resultMat) )
  expect_equal( ncol(res$matrix), ncol(resultMat) )
  expect_identical(res$matrix, resultMat)
  expect_equal( res$rowsKept, rownames(resultMat) )
  expect_equal( res$rowsRemoved, rownames(testMat)[3] )
  expect_equal( res$colsKept, colnames(resultMat) )
  expect_equal( res$colsRemoved, NULL )
})

test_that("check returned matrix: Rows and Cols", {
  resultMat <- testMat[ c(1:2,4:5), c(1:5,7,9,10) ]
  res <- RemoveZeroVariance(testMat, cols=TRUE)
  expect_equal( nrow(res$matrix), nrow(resultMat) )
  expect_equal( ncol(res$matrix), ncol(resultMat) )
  expect_identical(res$matrix, resultMat)
  expect_equal( res$rowsKept, rownames(resultMat) )
  expect_equal( res$rowsRemoved, rownames(testMat)[3] )
  expect_equal( res$colsKept, colnames(resultMat) )
  expect_equal( res$colsRemoved, colnames(testMat)[c(6,8)] )
})

test_that("check returned matrix: Cols Only", {
  resultMat <- testMat[ , c(1:5,7,9,10) ]
  res <- RemoveZeroVariance(testMat, rows=FALSE, cols=TRUE)
  expect_equal( nrow(res$matrix), nrow(resultMat) )
  expect_equal( ncol(res$matrix), ncol(resultMat) )
  expect_identical(res$matrix, resultMat)
  expect_equal( res$rowsKept, rownames(resultMat) )
  expect_equal( res$rowsRemoved, NULL )
  expect_equal( res$colsKept, colnames(resultMat) )
  expect_equal( res$colsRemoved, colnames(testMat)[c(6,8)] )
})


