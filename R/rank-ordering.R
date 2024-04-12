#' Rank order row data
#' 
#' Given a data-matrix of numeric data, calculates the rank of each row in each
#' column (feature in sample), gets the median rank across all columns, and
#' returns the original data with missing values set to NA, the reordered data,
#' and a data.frame of the ranks of each feature and the number of missing values.
#'
#' @param data_matrix matrix or data.frame of values
#' @param global_na the values to consider as missing
#' 
#' @export
#' 
#' @returns list with two matrices and a data.frame
rank_order_data = function(data_matrix, global_na = c(NA, Inf, 0))
{
  if (inherits(data_matrix, "data.frame")) {
    data_matrix = as.matrix(data_matrix)
  }
  missing_loc = setup_missing_matrix(data_matrix, global_na)
  data_matrix_na = data_matrix
  data_matrix_na[missing_loc] = NA
  sample_ranks = purrr::map(seq_len(ncol(data_matrix_na)), \(in_col){
    rank(data_matrix_na[, in_col], na.last = FALSE)
    
  })
  sample_ranks = do.call(cbind, sample_ranks)
  median_rank = apply(sample_ranks, 1, median)
  
  n_na = rowSums(is.na(data_matrix_na))
  rank_order = order(median_rank, decreasing = TRUE)
  
  perc_missing = colSums(is.na(data_matrix_na)) / nrow(data_matrix_na)
  perc_order = order(perc_missing, decreasing = TRUE)
  
  return(list(original = data_matrix_na,
              ordered = data_matrix_na[rank_order, perc_order],
              n_na_rank = data.frame(n_na = n_na,
                                     median_rank = median_rank)))
}
