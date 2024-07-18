#' Information-content-informed kendall tau
#' 
#' Given a data-matrix, computes the information-content-informed (ICI) Kendall-tau-b between
#' all samples.
#' 
#' @param data_matrix samples are columns, features are rows
#' @param global_na what values should be treated as missing (NA)?
#' @param zero_value what is the actual zero value?
#' @param perspective how to treat missing data in denominator and ties, see details
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_good should the diagonal entries reflect how many entries in the sample were "good"?
#' @param progress should progress be displayed.
#' 
#' @details For more details, see the ICI-Kendall-tau vignette:
#'   \href{../doc/ici-kendalltau.html}{\code{vignette("ici-kendalltau", package = "ICIKendallTau")}}
#' 
#' @return numeric
#' @keywords internal
#' 
ici_kendalltau_ref = function(data_matrix, 
                              global_na = c(NA, Inf, 0),
                              zero_value = 0, 
                              perspective = "global",
                              scale_max = TRUE,
                              diag_good = TRUE,
                              progress = FALSE){
  
  
  exclude_loc = matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  
  # Actual NA and Inf values are special cases, so we do
  # this very specifically.
  if (length(global_na) > 0) {
    if (any(is.na(global_na))) {
      exclude_loc[is.na(data_matrix)] = TRUE
      global_na = global_na[!is.na(global_na)]
    }
    
    if (any(is.infinite(global_na))) {
      exclude_loc[is.infinite(data_matrix)] = TRUE
      global_na = global_na[!is.infinite(global_na)]
    }
    
  }
  
  # now that we've done the NA and Inf values, we can go
  # ahead and take care of the rest.
  if (length(global_na) > 0) {
    for (ival in global_na) {
      exclude_loc[data_matrix == ival] = TRUE
    }
  }
  
  
  exclude_data = data_matrix
  exclude_data[exclude_loc] = NA
  # set everything to NA and let R take care of it
  
  cor_matrix = matrix(NA, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
  rownames(cor_matrix) = colnames(cor_matrix) = colnames(exclude_data)
  ntotal = 0
  for (icol in seq(1, ncol(exclude_data))) {
    for (jcol in seq(icol, ncol(exclude_data))) {
      cor_matrix[icol, jcol] = cor_matrix[jcol, icol] = ici_kt_pairs(exclude_data[, icol], exclude_data[, jcol], perspective = perspective)
      # ntotal = ntotal + 1
      # message(ntotal)
    }
  }
  
  # calculate the max-cor value for use in scaling across multiple comparisons
  n_observations = nrow(exclude_data)
  n_na = sort(colSums(exclude_loc))
  m_value = floor(sum(n_na[1:2]) / 2)
  n_m = n_observations - m_value
  max_cor_denominator = choose(n_m, 2) + n_observations * m_value
  max_cor_numerator = choose(n_m, 2) + n_observations * m_value + choose(m_value, 2)
  max_cor = max_cor_denominator / max_cor_numerator
  
  if (scale_max) {
    out_matrix = cor_matrix / max_cor
  } else {
    out_matrix = cor_matrix
  }
  
  if (diag_good) {
    n_good = colSums(!exclude_loc)
    diag(out_matrix) = n_good / max(n_good)
  }
  
  return(list(cor = out_matrix, raw = cor_matrix, keep = t(!exclude_loc)))
}

