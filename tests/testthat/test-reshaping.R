test_that("reshaping works properly", {
  n_val = 1000
  ind_elements = paste0("s", seq_len(n_val))
  cor_matrix_in = matrix(rnorm(n_val * n_val), nrow = n_val)
  rownames(cor_matrix_in) = colnames(cor_matrix_in) = ind_elements
  
  
  # set up a full one going from matrix to DF
  cor_df = cor_matrix_2_long_df(cor_matrix_in)
  expect_equal(nrow(cor_df), nrow(cor_matrix_in) * ncol(cor_matrix_in))
  expect_equal(cor_matrix_in["s2", "s342"], cor_df |> dplyr::filter(s1 %in% "s2", s2 %in% "s342") |> dplyr::pull(cor))
  
  # and then back again
  cor_matrix_out = long_df_2_cor_matrix(cor_df)
  expect_equal(cor_matrix_in["s2", "s342"], cor_matrix_out["s2", "s342"])
  
  # now setup a df that has only 1/2 the comparisons
  comb_elements = combn(ind_elements, 2)
  cor_df_short = data.frame(s1 = comb_elements[1, ],
                            s2 = comb_elements[2, ],
                            cor = cor_matrix_in[t(comb_elements)])
  cor_matrix_short = long_df_2_cor_matrix(cor_df_short)
  expect_equal(cor_matrix_in["s2", "s342"], cor_matrix_short["s2", "s342"])
  
  # and now don't expect a square output
  cor_matrix_ns = long_df_2_cor_matrix(cor_df_short, is_square = FALSE)
  expect_equal(nrow(cor_matrix_ns), n_val - 1)
  expect_equal(ncol(cor_matrix_ns), n_val - 1)
  expect_equal(cor_matrix_in["s2", "s342"], cor_matrix_ns["s2", "s342"])
  
  # check that we get an error if we don't have the right names
  cor_df_error = cor_df_short
  cor_df_error$raw = cor_df_error$cor
  expect_error(long_df_2_cor_matrix(cor_df_error[, c("s1", "s2", "raw")]), "must contain the names")
 })
