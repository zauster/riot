context("Testing GRAS algorithm")

A <- matrix(c(7, 3, 5, -3, 2, 9, 8, 0, -2, 0, 2, 0),
            ncol = 4, nrow = 3, byrow = TRUE)
u <- c(15, 26, -1)
v <- c(9, 16, 17, -2)

res <- doGRAS(A, u, v)


test_that("output is a matrix", {
    expect_true(is.matrix(res))
})

test_that("output matches", {

    ## first row
    expect_equal( res[1, 1], 8.424340881, tolerance = .002 )
  expect_equal( res[1, 2], 3.375237362, tolerance = .002 )
  expect_equal( res[1, 3], 5.200422021, tolerance = .002 )
  expect_equal( res[1, 4], -2.000000263, tolerance = .002 )

  ## second row
  expect_equal( res[2, 1], 3.000997078, tolerance = .002 )
  expect_equal( res[2, 2], 12.624763744, tolerance = .002 )
  expect_equal( res[2, 3], 10.374239178, tolerance = .002 )
  expect_equal( res[2, 4], 0, tolerance = .002 )

  ## third row
  expect_equal( res[3, 1],  -2.425339135, tolerance = .002 )
  expect_equal( res[3, 2], 0, tolerance = .002 )
  expect_equal( res[3, 3], 1.425339135, tolerance = .002 )
  expect_equal( res[3, 4], 0, tolerance = .002 )
})
