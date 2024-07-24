run_long_kendallt = as.logical(Sys.getenv("run_long_kendallt"))
if (is.na(run_long_kendallt)) {
  run_long_kendallt = FALSE
}
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
  
  expect_equal(ici_kt(x, y, perspective = "global")[c(1, 2)], ICIKendallTau:::ici_kt_pairs(x, y, "global"))
})

test_that("bad values passed error or return NA", {
  x = sort(rnorm(100))
  y = rep(NA, 100)
  result = numeric(3)
  result[1] = NA
  result[2] = NA
  result[3] = NA
  names(result) = c("tau", "pvalue", "tau_max")
  expect_equal(ici_kt(x, y), result)
  
  y = x[1:99]
  expect_error(ici_kt(x, y), "not the same length")
  
  expect_warning(ici_kt(x[2], y[2]), "vectors only have a single value")
  y2 = rep(1, 100)
  expect_warning(ici_kt(x, y2), "have only a single unique value")
})

test_that("matrix kendall works", {
  x = sort(rnorm(100))
  y = x + 1
  y[1:20] = NA
  
  test_mat = cbind(x, y)
  matrix_cor = ici_kendalltau(test_mat, global_na = c(NA),
                                  perspective = "global", scale_max = FALSE)
  expect_equal(ici_kt(x, y, "global")[[1]], matrix_cor$raw[2, 1])
})

test_that("large kendall returns correct", {
  set.seed(1234)
  x = rnorm(50000)
  y = rnorm(50000)
  ici_val = ici_kt(x, y, perspective = "global")
  expect_snapshot(ici_val)
})

if (run_long_kendallt) {
  test_that("big kendall works", {
    x = sort(rnorm(50000))
    y = x + 1
    x[1:5000] = NA
    
    t1 = ici_kt(x, y, perspective = "global")
    t2 = ICIKendallTau:::ici_kt_pairs(x, y, perspective = "global")
    expect_equal(t1[c(1, 2)], t2)
  })

  test_that("lots of kendall works", {
    x = matrix(rnorm(500000), 5000, 100)
    colnames(x) = paste0('s', seq_len(ncol(x)))

    t1 = ici_kendalltau(x)
    expect_lt(t1$run_time, 9)
  })
}


test_that("include_only works as intended", {
  set.seed(1234)
  x = matrix(rnorm(5000), nrow = 50, ncol = 100)
  colnames(x) = paste0("s", seq(1, ncol(x)))
  
  include_test = "s1"
  small_include = ici_kendalltau(x, include_only = include_test)
  expect_equal(sum(small_include$cor == 0), 9702)
  
  include_test = c("s1", "s3")
  small_include2 = ici_kendalltau(x, include_only = include_test)
  expect_equal(sum(small_include2$cor == 0), 9506)
  
  include_test = list(s1 = "s1", s2 = c("s2", "s3"))
  small_include3 = ici_kendalltau(x, include_only = include_test)
  expect_equal(sum(small_include3$cor == 0), 9896)
  
  # transform to a data.frame and make sure it runs the same
  include_df = as.data.frame(include_test)
  small_include4 = ici_kendalltau(x, include_only = include_df)
  expect_equal(small_include4$cor, small_include3$cor)
  
  # turn off self-comparisons and see if it gets smaller
  small_include5 = ici_kendalltau(x, include_only = include_test, diag_good = FALSE)
  expect_equal(sum(small_include5$cor == 0), 9996)
  
  small_include6 = ici_kendalltau(x, include_only = include_test, diag_good = FALSE, return_matrix = FALSE)
  expect_equal(nrow(small_include6$cor), 2)
})

test_that("completeness works correctly",{
  set.seed(1234)
  x = matrix(rnorm(5000), nrow = 50, ncol = 100, byrow = TRUE)
  colnames(x) = paste0("s", seq(1, ncol(x)))
  
  x[sample(5000, 40)] = NA
  x_cor = ici_kendalltau(x, perspective = "global", return_matrix = FALSE)
  x_comp = pairwise_completeness(x, return_matrix = FALSE)
  expect_equal(nrow(x_cor$cor), nrow(x_comp))
  
  expect_snapshot(x_comp[4:6, ])
})

test_that("kt_fast works properly", {
  set.seed(1234)
  x = matrix(rnorm(400), nrow = 100, ncol = 4)
  colnames(x) = paste0("s", seq(1, ncol(x)))
  
  fast_vals = kt_fast(x)
  base_vals = cor(x, method = "kendall")
  expect_equal(fast_vals$tau, base_vals)
  
  x_na = x
  x_na[, 1] = NA
  na_complete_1 = kt_fast(x_na, use = "complete.obs")
  expect_true(all(is.na(na_complete_1$tau)))
  
  x_na2 = x
  x_na2[10, 1] = NA
  
  na_pairs_everything = kt_fast(x_na2[, 1], x_na2[, 2])
  expect_equal(na_pairs_everything$tau[1, ], c('x_na2...1.' = NA_real_, 'x_na2...2.' = NA_real_))
  expect_snapshot(na_pairs_everything[c("tau", "pvalue")])
  na_pairs_complete = kt_fast(x_na2[, 1], x_na2[, 2], use = "complete.obs")
  expect_snapshot(na_pairs_complete[c("tau", "pvalue")])
  expect_gt(fast_vals$tau[1, 2], na_pairs_complete$tau[1, 2])
  
  na_pairs_pairwise = kt_fast(x_na2[, 1], x_na2[, 2], use = "pairwise.complete.obs")
  expect_equal(na_pairs_pairwise[c("tau", "pvalue")], na_pairs_complete[c("tau", "pvalue")])
  
  na_matrix_everything = kt_fast(x_na2, use = "everything")
  expect_equal(na_matrix_everything$tau[1, 1], NA_real_)
  expect_equal(na_matrix_everything$pvalue[1, 1], NA_real_)
  
  na_matrix_complete = kt_fast(x_na2, use = "complete.obs")
  expect_snapshot(na_matrix_complete[c("tau", "pvalue")])
  expect_lt(na_matrix_complete$tau[1, 2], fast_vals$tau[1, 2])
  
  na_matrix_pairwise = kt_fast(x_na2, use = "pairwise.complete.obs")
  expect_snapshot(na_matrix_pairwise[c("tau", "pvalue")])
  expect_equal(na_matrix_pairwise$tau[, 1], na_matrix_complete$tau[, 1])
  expect_equal(na_matrix_pairwise$tau[2, 1], na_pairs_complete$tau[2, 1])
  
  df_vals = kt_fast(as.data.frame(x))
  expect_equal(df_vals[c("tau", "pvalue")], fast_vals[c("tau", "pvalue")])
  
  df_out = kt_fast(x, return_matrix = FALSE)
  expect_equal(df_out$tau[4, "tau"], fast_vals$tau["s2", "s3"])
})

test_that("kt_fast warnings and errors come up", {
  set.seed(1234)
  x = matrix(rnorm(400), nrow = 100, ncol = 4)
  colnames(x) = paste0("s", seq(1, ncol(x)))
  
  expect_error(kt_fast(x, use = "na.or.complete"), "is not a supported")
  
  colnames(x) = NULL
  expect_error(kt_fast(x), 'Colnames of `x` must be be specified.')
  x_vec = x[, 1]
  expect_error(kt_fast(x_vec), '`x_vec` and `y` should both be provided as vectors, or `x_vec` should be matrix-like.')

  colnames(x) = paste0("s", seq(1, ncol(x)))
  y = x
  expect_error(kt_fast(x, y), 'Both `x` and `y` must be vectors.')

  x_vec = x[, 1]
  y_vec = y[, 2]

  tmp = kt_fast(x_vec, y_vec)
  expect_equal(colnames(tmp$tau), c("x_vec", "y_vec"))
})

test_that("big data kendall works", {
  x = sort(rnorm(100))
  y = x + 1
  y[1:20] = NA
  
  test_mat = cbind(x, y)
  matrix_cor = ici_kendalltau(test_mat, global_na = c(NA),
                              perspective = "global", scale_max = FALSE)
  long_cor = ici_kendalltau(test_mat, global_na = c(NA),
                              perspective = "global", scale_max = FALSE,
                              return_matrix = FALSE)
  
  expect_equal(nrow(long_cor$cor), 3)
  expect_equal(long_cor$cor$raw[1], matrix_cor$raw[2, 1])
  expect_equal(long_cor$cor$raw[3], matrix_cor$raw[2, 2])
  expect_equal(ici_kt(x, y, "global")[[1]], matrix_cor$raw[2, 1])
})

test_that("errors and messages appear", {
  x = matrix(rnorm(200), 20, 10)
  expect_error(ici_kendalltau(x), 'Colnames of `x` must be be specified.')
  expect_error(pairwise_completeness(x), 'Colnames of `x` must be be specified.')

  x_char = matrix(as.character(x), 20, 10)
  colnames(x_char) = paste0("S", seq_len(ncol(x)))
  expect_error(ici_kendalltau(x_char), '`x_char` must be a numeric type.')
  expect_error(pairwise_completeness(x_char), '`x_char` must be a numeric type.')

  x_df = as.data.frame(x)
  expect_message(ici_kendalltau(x_df), '`x_df` is a data.frame, converting to matrix ...')
  expect_message(pairwise_completeness(x_df), '`x_df` is a data.frame, converting to matrix ...')
})
