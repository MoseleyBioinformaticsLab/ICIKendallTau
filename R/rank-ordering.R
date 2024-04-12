#' Rank order row data
#' 
#' Given a data-matrix of numeric data, calculates the rank of each row in each
#' column (feature in sample), gets the median rank across all columns, and
#' returns the original data with missing values set to NA, the reordered data,
#' and a data.frame of the ranks of each feature and the number of missing values.
#'
#' @param data_matrix matrix or data.frame of values
#' @param global_na the values to consider as missing
#' @param sample_classes are the columns defined by some metadata?
#' 
#' @export
#' 
#' @returns list with two matrices and a data.frame
rank_order_data = function(data_matrix, global_na = c(NA, Inf, 0), 
                           sample_classes = NULL)
{
  if (inherits(data_matrix, "data.frame")) {
    data_matrix = as.matrix(data_matrix)
  }
  missing_loc = setup_missing_matrix(data_matrix, global_na)
  data_matrix_na = data_matrix
  data_matrix_na[missing_loc] = NA
  
  if (is.null(sample_classes)) {
    sample_classes = rep("rmf_abcd", ncol(data_matrix_na))
  }
  
  split_classes = split(colnames(data_matrix_na), sample_classes)
  
  get_ranks = function(in_na)
  {
    sample_ranks = purrr::map(seq_len(ncol(in_na)), \(in_col){
      rank(in_na[, in_col], na.last = FALSE)
    })
    sample_ranks = do.call(cbind, sample_ranks)
    
    median_rank = apply(sample_ranks, 1, median)
    n_na = rowSums(is.na(in_na))
    rank_order = order(median_rank, decreasing = TRUE)
    
    perc_missing = colSums(is.na(in_na)) / nrow(in_na)
    perc_order = order(perc_missing, decreasing = TRUE)
    list(original = in_na,
         ordered = in_na[rank_order, perc_order],
         n_na_rank = data.frame(n_na = n_na,
                                median_rank = median_rank))
  }
  
  split_ranks = purrr::imap(split_classes, \(in_split, split_id){
    split_na = data_matrix_na[, in_split]
    n_na = rowSums(is.na(split_na))
    keep_na = !(n_na == ncol(split_na))
    split_na = split_na[keep_na, ]
    
    na_info = get_ranks(split_na)
    if (!(split_id %in% "rmf_abcd")) {
      na_info$n_na_rank$split = split_id
    }
    na_info
  })
  
  if (length(split_ranks) == 1) {
    split_ranks = split_ranks[[1]]
  }
  
  return(split_ranks)
}
