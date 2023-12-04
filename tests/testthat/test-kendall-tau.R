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
  
  test_mat = rbind(x, y)
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
}

test_that("include_only works as intended", {
  set.seed(1234)
  x = matrix(rnorm(5000), nrow = 100, ncol = 50)
  rownames(x) = paste0("s", seq(1, nrow(x)))
  
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
})

test_that("completeness works correctly",{
  set.seed(1234)
  x = matrix(rnorm(5000), nrow = 100, ncol = 50)
  rownames(x) = paste0("s", seq(1, nrow(x)))
  
  x[sample(5000, 40)] = NA
  x_cor = ici_kendalltau(x, perspective = "global", return_matrix = FALSE)
  x_comp = pairwise_completeness(x, return_matrix = FALSE)
  expect_equal(nrow(x_cor$cor), nrow(x_comp))
  
  expect_snapshot(x_comp[4:6, ])
})

test_that("fast_kt works properly", {
  set.seed(1234)
  x = matrix(rnorm(400), nrow = 100, ncol = 4)
  colnames(x) = paste0("s", seq(1, ncol(x)))
  
  fast_vals = kt_fast(x)
  base_vals = cor(x, method = "kendall")
  expect_equal(fast_vals$tau, base_vals)
  
  expect_error(kt_fast(x, use = "na.or.complete"), "is not a supported")
  
  x_na = x
  x_na[, 1] = NA
  na_complete_1 = kt_fast(x_na, use = "complete.obs")
  expect_true(all(is.na(na_complete_1$tau)))
  
  x_na2 = x
  x_na2[10, 1] = NA
  
  na_pairs_everything = kt_fast(x_na2[, 1], x_na2[, 2])
  expect_equal(na_pairs_everything, c("tau" = NA, "pvalue" = NA))
  expect_snapshot(na_pairs_everything[c("tau", "pvalue")])
  na_pairs_complete = kt_fast(x_na2[, 1], x_na2[, 2], use = "complete.obs")
  expect_snapshot(na_pairs_complete[c("tau", "pvalue")])
  expect_gt(fast_vals$tau[1, 1], na_pairs_complete["tau"])
  
  na_pairs_pairwise = kt_fast(x_na2[, 1], x_na2[, 2], use = "pairwise.complete.obs")
  expect_equal(na_pairs_pairwise, na_pairs_complete)
  
  na_matrix_everything = kt_fast(x_na2, use = "everything")
  expect_equal(na_matrix_everything$tau[1, 1], as.double(NA))
  expect_equal(na_matrix_everything$pvalue[1, 1], as.double(NA))
  expect_equal(na_matrix_everything$tau[3, 2], fast_vals$tau[3, 2])
  expect_equal(na_matrix_everything$pvalue[3, 2], fast_vals$pvalue[3, 2])

  
  na_matrix_complete = kt_fast(x_na2, use = "complete.obs")
  expect_snapshot(na_matrix_complete[c("tau", "pvalue")])
  expect_lt(na_matrix_complete$tau[1, 2], fast_vals$tau[1, 2])
  
  na_matrix_pairwise = kt_fast(x_na2, use = "pairwise.complete.obs")
  expect_snapshot(na_matrix_pairwise[c("tau", "pvalue")])
  expect_equal(na_matrix_pairwise$tau[, 1], na_matrix_complete$tau[, 1])
  expect_equal(na_matrix_pairwise$tau[2, 3], na_matrix_everything$tau[2, 3])
  expect_equal(na_matrix_pairwise$tau[2, 1], na_pairs_complete[["tau"]])
  
  df_vals = kt_fast(as.data.frame(x))
  expect_equal(df_vals[c("tau", "pvalue")], fast_vals[c("tau", "pvalue")])
  
  df_out = kt_fast(x, return_matrix = FALSE)
  expect_equal(df_out$tau[4, "tau"], fast_vals$tau["s2", "s3"])
})
