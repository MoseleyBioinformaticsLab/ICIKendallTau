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

test_that("matrix kendall works", {
  x = sort(rnorm(100))
  y = x + 1
  y[1:20] = NA
  
  test_mat = rbind(x, y)
  matrix_cor = ici_kendalltau(test_mat, exclude_na = TRUE, exclude_0 = FALSE,
                                  perspective = "global", scale_max = FALSE)
  expect_equal(ici_kt(x, y, "global")[[1]], matrix_cor$cor[2, 1])
})
