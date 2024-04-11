rank_order_data = function(data_matrix, global_na = c(NA, Inf, 0))
{
  missing_loc = setup_missing_matrix(data_matrix, global_na)
  data_matrix_na = data_matrix
  data_matrix_na[missing_loc] = NA
  sample_ranks = purrr::map(seq_len(ncol(data_matrix_na)), \(in_col){
    rank(data_matrix_na[, in_col], na.last = FALSE)
    
  })
  sample_ranks = do.call(cbind, sample_ranks)
  median_ranks = apply(sample_ranks, 1, median)
  rank_order = order(median_ranks, decreasing = TRUE)
  
  perc_missing = colSums(is.na(data_matrix_na)) / nrow(data_matrix_na)
  perc_order = order(perc_missing, decreasing = TRUE)
  
  return(original = data_matrix_na,
         ordered = data_matrix_na[rank_order, perc_order])
}
