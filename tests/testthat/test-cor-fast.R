test_that("other correlations work", {
  set.seed(1234)
  x = matrix(rnorm(400), nrow = 100, ncol = 4)
  colnames(x) = paste0("s", seq(1, ncol(x)))

  p_res = cor_fast(x, method = "pearson")
  s_res = cor_fast(x, method = "spearman")

  expect_gt(s_res$pvalue[2, 1], p_res$pvalue[2, 1])

  expect_equal(ncol(p_res$pvalue), 4)

  base_pearson = cor(x, method = "pearson")
  test_pearson = stats::cor.test(x[, 1], x[, 2], method = "pearson")
  
  expect_equal(p_res$rho, base_pearson)
  expect_equal(p_res$rho[2, 1], test_pearson[["estimate"]][[1]])
  expect_equal(p_res$pvalue[2, 1], test_pearson[["p.value"]])
  
  base_spearman = cor(x, method = "spearman")
  test_spearman = stats::cor.test(x[, 1], x[, 2], method = "spearman")
  expect_equal(s_res$rho, base_spearman)
  expect_equal(s_res$pvalue[2, 1], test_spearman[["p.value"]])
})

test_that("cor_fast warnings and errors come up", {
  set.seed(1234)
  x = matrix(rnorm(400), nrow = 100, ncol = 4)
  colnames(x) = paste0("s", seq(1, ncol(x)))
  
  expect_error(cor_fast(x, use = "na.or.complete"), "is not a supported")
  
  colnames(x) = NULL
  expect_error(cor_fast(x), 'Colnames of `x` must be be specified.')
  x_vec = x[, 1]
  expect_error(cor_fast(x_vec), '`x_vec` and `y` should both be provided as vectors, or `x_vec` should be matrix-like.')

  colnames(x) = paste0("s", seq(1, ncol(x)))
  y = x
  expect_error(cor_fast(x, y), 'Both `x` and `y` must be vectors.')

  x_vec = x[, 1]
  y_vec = y[, 2]

  tmp = cor_fast(x_vec, y_vec)
  expect_equal(colnames(tmp$rho), c("x_vec", "y_vec"))
})

test_that("na values are handled correctly by cor_fast", {
  set.seed(1234)
  x = matrix(rnorm(400), nrow = 100, ncol = 4)
  colnames(x) = paste0("s", seq(1, ncol(x)))

  na_locs = sample(400, 40)
  x[na_locs] = NA

  cor_base_complete = cor.test(x[, 1], x[, 2], method = "pearson", use = "complete")
  cor_fast_complete = cor_fast(x = x[, 1], y = x[, 2], method = "pearson", use = "complete", return_matrix = FALSE)
  expect_equal(cor_base_complete$estimate[["cor"]], cor_fast_complete$rho$rho[1])

  cor_fast_pairwise_1 = cor_fast(x = x[, 1], y = x[, 2], method = "pearson", use = "pairwise.complete.obs", return_matrix = FALSE)
  cor_fast_pairwise_all = cor_fast(x, method = "pearson", use = "pairwise.complete.obs", return_matrix = FALSE)
  expect_equal(cor_fast_pairwise_1$rho$rho[1], cor_fast_pairwise_all$rho$rho[1])
  
  cor_base_pairwise_2 = cor.test(x[, 1], x[, 3], method = "pearson", use = "pairwise.complete.obs")
  expect_equal(cor_fast_pairwise_all$rho$rho[2], cor_base_pairwise_2$estimate[["cor"]])
})