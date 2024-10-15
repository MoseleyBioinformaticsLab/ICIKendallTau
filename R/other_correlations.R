#' Fast correlation with test
#' 
#' Allows to run `cor.test` on a matrix of inputs.
#' 
#' @param x a numeric vector, matrix, or data frame.
#' @param y NULL (default) or a vector.
#' @param use an optional character string giving a method for computing correlations in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", or "pairwise.complete.obs".
#' @param method which correlation method to use, "pearson" or "spearman"
#' @param alternative how to perform the statistical test
#' @param continuity should a continuity correction be applied
#' @param return_matrix should the matrices of values be returned, or a long data.frame
#' 
#' @details Although the interface is *mostly* identical to the built-in 
#' [stats::cor.test()] method, there are some differences. 
#'   
#'  * if only `x` is provided as a matrix, the columns must be named.
#'  * if providing both `x` and `y`, it is assumed they are both
#'  single vectors.
#'  * if `NA` values are present, this function does not error, but will either remove them 
#'  or return `NA`, depending on the option.
#'  * "na.or.complete" is not a valid option for `use`.
#'  * A named list with matrices or data.frame is returned, with the `rho` and `pvalue` values.
#' 
#' @return a list of matrices, rho, pvalue, or a data.frame.
#' @export
cor_fast = function(x, y = NULL, use = "everything", 
                    method = "pearson", 
                    alternative = "two.sided", 
                    continuity = FALSE,
                    return_matrix = TRUE)
{
  do_log_memory = get("memory", envir = icikt_logger)
  use_arg = rlang::caller_arg(use)
  x_arg = rlang::caller_arg(x)
  y_arg = rlang::caller_arg(y)
  
  # checking use
  na_method = match.arg(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                               "everything", "na.or.complete"))
  
  #message(na_method)
  # all of this top stuff can probably be wrapped in a function to parse
  # the logic, and can borrow some of the other new helpers.
  check_na_method(na_method)

  x = check_x_y(x = x, y = y,
            x_arg = x_arg, 
            y_arg = y_arg)

  na_vals = is.na(x)
  any_na = any(na_vals)
  no_na_rows = rowSums(na_vals) == 0
  sum_nona = sum(no_na_rows)

  computation = check_furrr()

  split_comparisons = setup_comparisons(colnames(x),
                                        include_only = NULL,
                                        diag_good = FALSE,
                                        ncore = computation$ncore,
                                        which = method,
                                        include_arg = NULL)
  
  do_computation = TRUE
  if (na_method %in% c("everything", "all.obs")) {
    if (any_na) {
      do_computation = FALSE
    }
  }

  if (na_method %in% c("complete.obs", "pariwise.complete.obs")) {
    if (sum_nona == 0) {
      do_computation = FALSE
    } else {
      x = x[no_na_rows, , drop = FALSE]
    }
  }

  if (!do_computation) {
    split_cor = purrr::map(split_comparisons, \(in_split){
      in_split$rho = NA
      in_split$pvalue = NA
      in_split
    })
  }

  
  
  # we record how much time is actually spent doing ICI-Kt
  # itself, as some of the other operations will add a bit of time
  # 
  if (do_computation) {
    t1 = Sys.time()
    # note here, this takes our list of comparisons, and then calls the do_split
    # function above on each of them.
    split_cor = computation$split_fun(split_comparisons, cor_split, x, na_method, do_log_memory, method, alternative, continuity)
    t2 = Sys.time()
    t_diff = as.numeric(difftime(t2, t1, units = "secs"))
  } else {
    t_diff = 0
  }
  
  all_cor = purrr::list_rbind(split_cor)
  if (return_matrix) {
    log_message("Generating the output matrix ...")
    cor_matrix = matrix(0, nrow = ncol(x), ncol = ncol(x))
    rownames(cor_matrix) = colnames(cor_matrix) = colnames(x)
    rho_matrix = cor_matrix
    pvalue_matrix = cor_matrix
    
    one_way_index = cbind(all_cor$s1, all_cor$s2)
    back_way_index = cbind(all_cor$s2, all_cor$s1)
    
    rho_matrix[one_way_index] = all_cor$rho
    rho_matrix[back_way_index] = all_cor$rho
    
    pvalue_matrix[one_way_index] = all_cor$pvalue
    pvalue_matrix[back_way_index] = all_cor$pvalue
    
    return(list(rho = rho_matrix, pvalue = pvalue_matrix, run_time = t_diff))
  } else {
    return(list(rho = all_cor, run_time = t_diff))
  }
  
}

cor_split = function(do_comparisons, x, na_method, do_log_memory, method, alternative, continuity) {
  #seq_range = seq(in_range[1], in_range[2])
  #print(seq_range)
  rho = vector("numeric", nrow(do_comparisons))
  pvalue = rho
  
  for (irow in seq(1, nrow(do_comparisons))) {
    return_na = FALSE
    iloc = do_comparisons[irow, 1]
    jloc = do_comparisons[irow, 2]
    
    tmp_x = x[, iloc]
    tmp_y = x[, jloc]
    if (na_method %in% "pairwise.complete.obs") {
      pair_good = !is.na(tmp_x) & !is.na(tmp_y)
      if (sum(pair_good) == 0) {
        return_na = TRUE
      } else {
        return_na = FALSE
        tmp_x = tmp_x[pair_good]
        tmp_y = tmp_y[pair_good]
      }
    } else if (na_method %in% c("everything", "all.obs")) {
      return_na = any(is.na(c(tmp_x, tmp_y)))
    }
    
    if (return_na) {
      rho[irow] = as.double(NA)
      pvalue[irow] = as.double(NA)
      
    } else {
      cor_res = stats::cor.test(tmp_x, tmp_y,
                                alternative = alternative,
                                method = method,
                                continuity = continuity)
      rho[irow] = cor_res[["estimate"]][[1]]
      pvalue[irow] = cor_res[["p.value"]]
      
    }
    if (do_log_memory && ((irow %% 100) == 0)) {
      log_memory()
    }
    
  }
  do_comparisons$rho = rho
  do_comparisons$pvalue = pvalue
  do_comparisons
}
