test_that("basic kendall-tau matches base R", {
  x = seq(1, 10)
  y = seq(1, 10)
  expect_equal(ici_kt(x, y)[[1]], cor(x, y, method = "kendall"))
  
  y[2] = 15
  expect_equal(ici_kt(x, y)[[1]], cor(x, y, method = "kendall"))
  
  y = seq(10, 1) # should give -1
  expect_equal(ici_kt(x, y)[[1]], cor(x, y, method = "kendall"))
  
  y[2] = 15
  expect_equal(ici_kt(x, y)[[1]], cor(x, y, method = "kendall"))
})

test_that("difference and reference match - short", {
  x = sort(rnorm(100))
  y = x + 1
  y[1:20] = NA
  
  expect_equal(ici_kt(x, y, "global"), ICIKendallTau:::ici_kt_pairs(x, y, "global"))
})

test_that("bad values passed return NA", {
  x = sort(rnorm(100))
  y = rep(NA, 100)
  result = numeric(2)
  result[1] = NA
  result[2] = NA
  names(result) = c("tau", "pvalue")
  expect_equal(ici_kt(x, y), result)
  
  y = x[1:99]
  expect_error(ici_kt(x, y), "not the same length")
  
  expect_equal(ici_kt(x[1], y[1]), result)
})

test_that("matrix kendall works", {
  x = sort(rnorm(100))
  y = x + 1
  y[1:20] = NA
  
  test_mat = rbind(x, y)
  matrix_cor = ici_kendalltau(test_mat, exclude_na = TRUE, exclude_0 = FALSE,
                                  perspective = "global", scale_max = FALSE)
  expect_equal(ici_kt(x, y, "global")[[1]], matrix_cor$cor[2, 1])
})
