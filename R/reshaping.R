utils::globalVariables(c("stack", "ind", "values"))

#' convert matrix to data.frame
#' 
#' Given a square correlation matrix, converts it to a long data.frame, with three columns.
#' 
#' @param in_matrix the correlation matrix
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
  
  wide_df = as.data.frame(in_matrix)
  wide_df$row = rownames(in_matrix)
  
  long_df = data.frame(row = wide_df$row, stack(wide_df, select = -row))
  long_df$s1 = as.character(long_df$row)
  long_df$s2 = as.character(long_df$ind)
  long_df$cor = long_df$values
  long_df$row = NULL
  long_df$ind = NULL
  long_df$values = NULL
  
  long_df
}

#' convert data.frame to matrix
#' 
#' Given a long data.frame, converts it to a possibly square correlation matrix
#' 
#' @param long_df the long data.frame
#' @param is_square should it be a square matrix?
#' 
#' 
#' @export
#' @return matrix
long_df_2_cor_matrix = function(long_df, is_square = TRUE)
{
  check_names = all(c("s1", "s2", "cor") %in% names(long_df))
  if (!check_names) {
    stop("The data.frame must contain the names 's1', 's2', and 'cor'.")
  }
  
  if (is_square) {
    s1_factors = s2_factors = factor(c(long_df[["s1"]], long_df[["s2"]]))
  } else {
    s1_factors = factor(long_df[["s1"]])
    s2_factors = factor(long_df[["s2"]])
  }
  cor_matrix = matrix(nrow = nlevels(s1_factors), ncol = nlevels(s2_factors),
                      dimnames = list(levels(s1_factors), levels(s2_factors)))
  cor_matrix[cbind(long_df[["s1"]], long_df[["s2"]])] = long_df[["cor"]]
  
  # if the number of rows aren't equal, then we only have one half of
  # the matrix, and we need to grab the rest.
  # Otherwise, we assume we have all of the matrix, and we are good.
  if ((nrow(long_df) != (nrow(cor_matrix) * ncol(cor_matrix))) && (is_square)) {
    cor_matrix[cbind(long_df[["s2"]], long_df[["s1"]])] = long_df[["cor"]]
  }
  
  return(cor_matrix)
}
