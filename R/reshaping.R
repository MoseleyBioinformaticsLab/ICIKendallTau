#' convert matrix to data.frame
#' 
#' Given a square correlation matrix, converts it to a long data.frame, with three columns.
#' 
#' @param cor_matrix the correlation matrix
#' 
#' @details The data.frame contains three columns:
#' * s1: the first entry of comparison
#' * s2: the second entry of comparison
#' * cor: the correlation value
#' 
#' @export
#' @return data.frame
cor_matrix_2_long_df = function(in_matrix)
{
  # in_matrix = pcor_vals
  if (nrow(in_matrix) != ncol(in_matrix)) {
    stop("Matrix isn't square, can't do it!")
  } 
  
  wide_df = as.data.frame(in_matrix)
  wide_df$row = rownames(in_matrix)
  
  if ("dplyr" %in% utils::installed.packages()) {
    long_df = data.frame(row = wide_df$row, stack(wide_df, select = -row)) |>
      dplyr::transmute(s1 = row, s2 = ind, cor = values)
  } else {
    stop("This operation uses 'dplyr', and that package is not installed.")
  }
  
  long_df
}

#' convert data.frame to matrix
#' 
#' Given a long data.frame, converts it to a square correlation matrix
#' 
#' @param long_df the long data.frame
#' 
#' 
#' @export
#' @return matrix
long_df_2_cor_matrix = function(long_df)
{
  check_names = all(names(long_df) %in% c("s1", "s2", "cor"))
  if (!check_names) {
    stop("The data.frame must contain the names 's1', 's2', and 'cor'.")
  }
  long_df$s1 = as.factor(long_df$s1)
  long_df$s2 = as.factor(long_df$s2)
  cor_matrix = matrix(nrow = nlevels(long_df$s1), ncol = nlevels(long_df$s2),
                      dimnames = list(levels(long_df$s1), levels(long_df$s2)))
  cor_matrix[cbind(long_df$s1, long_df$s2)] = long_df$cor
  cor_matrix[cbind(long_df$s2, long_df$s1)] = long_df$cor
  return(cor_matrix)
}
