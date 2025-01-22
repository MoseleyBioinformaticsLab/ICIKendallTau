#' Information-content-informed kendall tau
#' 
#' Given a data-matrix, computes the information-theoretic Kendall-tau-b between
#' all samples.
#' 
#' @param data_matrix matrix or data.frame of values, samples are columns, features are rows
#' @param global_na numeric vector that defines globally, what should be treated as NA?
#' @param perspective how to treat missing data in denominator and ties, character
#' @param scale_max logical, should everything be scaled compared to the maximum correlation?
#' @param diag_good logical, should the diagonal entries reflect how many entries in the sample were "good"?
#' @param include_only only run the correlations that include the members (as a vector) or combinations (as a list or data.frame)
#' @param alternative what is the alternative for the p-value test?
#' @param continuity should a continuity correction be applied?
#' @param check_timing logical to determine should we try to estimate run time for full dataset? (default is FALSE)
#' @param return_matrix logical, should the data.frame or matrix result be returned?
#' 
#' @seealso [test_left_censorship()] [pairwise_completeness()] [kt_fast()]
#' 
#' @details For more details, see the vignette `vignette("ici-kendalltau", package = "ICIKendallTau"))`
#' 
#'   The default for \code{global_na} includes what values in the data to replace with NA for the Kendall-tau calculation. By default these are \code{global_na = c(NA, Inf, 0)}. If you want to replace something other than 0, for example, you might use \code{global_na = c(NA, Inf, -2)}, and all values of -2 will be replaced instead of 0.
#'   
#'   When \code{check_timing = TRUE}, 5 random pairwise comparisons will be run to generate timings on a single core, and then estimates of how long the full set will take are calculated. The data is returned as a data.frame, and will be on the low side, but it should provide you with a good idea of how long your data will take.
#'   
#'   Returned is a list containing matrices with:
#'   
#'   * cor: scaled correlations
#'   * raw: raw kendall-tau correlations
#'   * pvalue: p-values
#'   * taumax: the theoretical maximum kendall-tau value possible
#'   * completeness: how complete the two samples are (i.e. how many entries are not missing in either sample)
#'   
#'   Eventually, we plan to provide two more parameters for replacing values, \code{feature_na} for feature specific NA values and \code{sample_na} for sample specific NA values.
#'   
#'   If you want to know if the missing values in your data are possibly due to 
#'   left-censorship, we recommend testing that hypothesis with [test_left_censorship()]
#'   first.
#' 
#' @return list with cor, raw, pvalue, taumax, completeness
#' 
#' @examples
#' \dontrun{
#' # not run
#' set.seed(1234)
#' s1 = sort(rnorm(1000, mean = 100, sd = 10))
#' s2 = s1 + 10 
#' 
#' matrix_1 = cbind(s1, s2)
#' 
#' r_1 = ici_kendalltau(matrix_1)
#' r_1$cor
#' 
#' #    s1 s2
#' # s1  1  1
#' # s2  1  1
#' names(r_1)
#' # "cor", "raw", "pvalue", "taumax", "completeness", "keep", "run_time"
#' 
#' s3 = s1
#' s3[sample(100, 50)] = NA
#' 
#' s4 = s2
#' s4[sample(100, 50)] = NA
#' 
#' matrix_2 = cbind(s3, s4)
#' r_2 = ici_kendalltau(matrix_2)
#' r_2$cor
#' #           s3        s4
#' # s3 1.0000000 0.9944616
#' # s4 0.9944616 1.0000000
#' 
#' # using include_only
#' set.seed(1234)
#' x = t(matrix(rnorm(5000), nrow = 100, ncol = 50))
#' colnames(x) = paste0("s", seq(1, nrow(x)))
#' 
#' # only calculate correlations of other columns with "s1"
#' include_s1 = "s1"
#' s1_only = ici_kendalltau(x, include_only = include_s1)
#' 
#' # include s1 and s3 things both
#' include_s1s3 = c("s1", "s3")
#' s1s3_only = ici_kendalltau(x, include_only = include_s1s3)
#' 
#' # only specify certain pairs either as a list
#' include_pairs = list(g1 = "s1", g2 = c("s2", "s3"))
#' s1_other = ici_kendalltau(x, include_only = include_pairs)
#' 
#' # or a data.frame
#' include_df = as.data.frame(list(g1 = "s1", g2 = c("s2", "s3")))
#' s1_df = ici_kendalltau(x, include_only = include_df)
#' 
#' }
#' @export
#' 
ici_kendalltau = function(data_matrix, 
                          global_na = c(NA, Inf, 0),
                          perspective = "global",
                          scale_max = TRUE,
                          diag_good = TRUE,
                          include_only = NULL,
                          alternative = "two.sided",
                          continuity = FALSE,
                          check_timing = FALSE,
                          return_matrix = TRUE){
  
  arg = rlang::caller_arg(data_matrix)
  include_arg = rlang::caller_arg(include_only)
  diag_arg = rlang::caller_arg(diag_good)
  
  do_log_memory = get("memory", envir = icikt_logger)
  
  check_if_colnames_null(data_matrix, arg = arg)
  data_matrix = transform_to_matrix(data_matrix, arg = arg)
  check_if_numeric(data_matrix, arg = arg)

  log_message("Processing missing values ...\n")
  
  exclude_loc = setup_missing_matrix(data_matrix, global_na)
  exclude_data = data_matrix
  exclude_data[exclude_loc] = NA
  n_sample = ncol(exclude_data)
  n_feature = nrow(exclude_data)
  
  
  # figure out if we are using furrr to process these
  computation = check_furrr()
  
  # generate the array of comparisons, 2 x ...,
  # where each column is a comparison between two columns of data
  log_message("Figuring out comparisons to do ...")
  split_comparisons = setup_comparisons(colnames(data_matrix),
                                        include_only = include_only,
                                        diag_good = diag_good,
                                        ncore = computation$ncore,
                                        which = "icikt",
                                        include_arg = include_arg)
  
  n_todo = sum(purrr::map_int(split_comparisons, nrow))
  
  if (check_timing) {
    
    sample_compare = sample(nrow(split_comparisons[[1]]), 5)
    tmp_pairwise = split_comparisons[[1]][sample_compare, , drop = FALSE]
    
    run_tmp = check_icikt_timing(exclude_data, tmp_pairwise, perspective, n_todo, computation$ncore)
    return(run_tmp)
  }
  
  
  # we record how much time is actually spent doing ICI-Kt
  # itself, as some of the other operations will add a bit of time
  # 
  log_message("Running correlations ...")
  t1 = Sys.time()
  # note here, this takes our list of comparisons, and then calls the do_split
  # function above on each of them.
  split_cor = computation$split_fun(split_comparisons, ici_split, exclude_data, perspective, do_log_memory, alternative, continuity)
  t2 = Sys.time()
  t_diff = as.numeric(difftime(t2, t1, units = "secs"))
  
  # put all the results back together again into one data.frame
  log_message("Recombining results ...")
  
  n_good = colSums(!exclude_loc)
  names(n_good) = colnames(data_matrix)
  frac_complete = n_good / nrow(exclude_loc)
  scaled_output = scale_and_reshape(split_cor,
                                    n_good = n_good,
                                    scale_max = scale_max,
                                    diag_good = diag_good,
                                    frac_complete = frac_complete,
                                    return_matrix = return_matrix,
                                    exclude_loc = exclude_loc) 
     
  scaled_output[["run_time"]] = t_diff
  return(scaled_output)
    
}

setup_comparisons = function(samples,
                              include_only = include_only,
                              diag_good = diag_good,
                              ncore = ncore,
                              which = "icikt",
                            include_arg = rlang::caller_arg(include_only))
{
  n_sample = length(samples)
  pairwise_comparisons = utils::combn(n_sample, 2)
  
  if (!diag_good) {
    extra_comparisons = matrix(rep(seq(1, n_sample), each = 2), nrow = 2, ncol = n_sample, byrow = FALSE)
    pairwise_comparisons = cbind(pairwise_comparisons, extra_comparisons)
  }
  
  # create a data.frame of the comparisons by the names of the columns instead,
  # this enables indexed comparisons and named comparisons, because
  # we can have row / column names in R
  # This is now n_comparisons x 2
  named_comparisons = data.frame(s1 = samples[pairwise_comparisons[1, ]],
                                 s2 = samples[pairwise_comparisons[2, ]])
  
  if (!is.null(include_only)) {
    if (is.character(include_only) || is.numeric(include_only)) {
      n_include = length(include_only)
      #message("a vector!")
      # Check each of the comparison vectors against the include_only variable
      # This returns TRUE where they match
      # Use OR to make sure we return everything that should be returned
      s1_include = named_comparisons$s1 %in% include_only
      s2_include = named_comparisons$s2 %in% include_only
      named_comparisons = named_comparisons[(s1_include | s2_include), ]
    } else if (is.list(include_only)) {
      if (length(include_only) == 2) {
        n_include = length(include_only[[1]])
        #message("a list!")
        # In this case the include_only is a list of things, so we have to check both
        # of the sets against each of the lists. Again, this returns TRUE where
        # they match. Because we want the things where they
        # are both TRUE (assuming l1[1] goes with l2[1]), we use the AND at the end.
        l1_l2 = paste0(include_only[[1]], "-", include_only[[2]])
        l2_l1 = paste0(include_only[[2]], "-", include_only[[1]])
        
        all_include = c(l1_l2, l2_l1)
        all_compare = paste0(named_comparisons$s1, "-", named_comparisons$s2)
        
        keep_compare = all_compare %in% all_include
        named_comparisons = named_comparisons[keep_compare, ]
        
      } else {
        cli::cli_abort(message = c(
          '{.arg {include_arg}} must be a vector, a data.frame with two columns, or list of two vectors.',
          'x' = 'Currently, {.code {length({include_arg})}} returns \\
          {length(include_only)}'
        ))
      }
    }
  }
  
  n_todo = nrow(named_comparisons)
  if (n_todo == 0) {
    cli::cli_abort(message = c(
      'No comparisons to do.',
      'i' = 'Check the list of column names in {.arg {include_arg}} vs those in the samples.'
    ))
    
  }
  
  log_message("Splitting up across compute ...")
  n_each = ceiling(n_todo / ncore)
  
  which_core = rep(seq(1, ncore), each = n_each)
  which_core = which_core[seq_len(nrow(named_comparisons))]
  
  named_comparisons$core = which_core
  if (which %in% "icikt") {
    named_comparisons$raw = Inf
    named_comparisons$pvalue = Inf
    named_comparisons$taumax = Inf
  } else if (which %in% "completeness") {
    named_comparisons$missingness = Inf
    named_comparisons$completeness = Inf
  } else if (which %in% "kt") {
    named_comparisons$tau = Inf
    named_comparisons$pvalue = Inf
  } else if (which %in% "pearson") {
    named_comparisons$rho = Inf
    named_comparisons$pvalue = Inf
    named_comparisons$n_values = Inf
  } else if (which %in% "spearman") {
    named_comparisons$rho = Inf
    named_comparisons$pvalue = Inf
    named_comparisons$n_values = Inf
  }
    
  split_comparisons = split(named_comparisons, named_comparisons$core)
  return(split_comparisons)
}

ici_split = function(do_comparisons, exclude_data, perspective, do_log_memory, alternative, continuity) {
  #seq_range = seq(in_range[1], in_range[2])
  #print(seq_range)
  
  raw = vector("numeric", nrow(do_comparisons))
  pvalue = raw
  taumax = raw
  completeness = raw
  
  for (irow in seq_len(nrow(do_comparisons))) {
    iloc = do_comparisons[irow, 1]
    jloc = do_comparisons[irow, 2]
    ici_res = ici_kt(exclude_data[, iloc], exclude_data[, jloc], perspective = perspective,
                      alternative = alternative, continuity = continuity)
    raw[irow] = ici_res["tau"]
    pvalue[irow] = ici_res["pvalue"]
    taumax[irow] = ici_res["tau_max"]
    completeness[irow] = ici_res["completeness"]
    if (do_log_memory && ((irow %% 100) == 0)) {
      log_memory()
    }
  }
  do_comparisons$raw = raw
  do_comparisons$pvalue = pvalue
  do_comparisons$taumax = taumax
  do_comparisons$completeness = completeness
  #return(ls())
  do_comparisons
}

kt_split = function(do_comparisons, x, na_method, do_log_memory) {
  #seq_range = seq(in_range[1], in_range[2])
  #print(seq_range)
  tau = vector("numeric", nrow(do_comparisons))
  pvalue = tau
  
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
      tau[irow] = as.double(NA)
      pvalue[irow] = as.double(NA)
      
    } else {
      ici_res = ici_kt(tmp_x, tmp_y)
      tau[irow] = ici_res["tau"]
      pvalue[irow] = ici_res["pvalue"]
      
    }
    if (do_log_memory && ((irow %% 100) == 0)) {
      log_memory()
    }
    
  }
  do_comparisons$tau = tau
  do_comparisons$pvalue = pvalue
  do_comparisons
}


scale_and_reshape = function(split_cor,
                              n_good = n_good,
                              scale_max = scale_max,
                              diag_good = diag_good,
                              frac_complete = frac_complete,
                              return_matrix = return_matrix,
                            exclude_loc = exclude_loc) 
{
  all_cor = purrr::list_rbind(split_cor)
  rownames(all_cor) = NULL

  if (scale_max) {
    max_cor = max(all_cor$taumax, na.rm = TRUE)
    all_cor$cor = all_cor$raw / max_cor
  } else {
    all_cor$cor = all_cor$raw
  }
  
  if (diag_good) {
    extra_cor = data.frame(s1 = names(n_good),
                           s2 = names(n_good),
                           core = 0,
                           raw = n_good / max(n_good),
                           pvalue = 0,
                           taumax = 1,
                           completeness = frac_complete,
                           cor = n_good / max(n_good))
    rownames(extra_cor) = NULL
    all_cor = rbind(all_cor, extra_cor)
  }
  
  # if the user asks for the matrix back, we give the matrices, otherwise we
  # leave them in the data.frame
  if (return_matrix) {
    log_message("Generating the output matrix ...")
    cor_matrix = matrix(0, nrow = length(n_good), ncol = length(n_good))
    rownames(cor_matrix) = colnames(cor_matrix) = names(n_good)
    raw_matrix = cor_matrix
    pvalue_matrix = cor_matrix
    taumax_matrix = cor_matrix
    completeness_matrix = cor_matrix
    
    one_way_index = cbind(all_cor$s1, all_cor$s2)
    back_way_index = cbind(all_cor$s2, all_cor$s1)
    
    raw_matrix[one_way_index] = all_cor$raw
    raw_matrix[back_way_index] = all_cor$raw
    
    cor_matrix[one_way_index] = all_cor$cor
    cor_matrix[back_way_index] = all_cor$cor
    
    pvalue_matrix[one_way_index] = all_cor$pvalue
    pvalue_matrix[back_way_index] = all_cor$pvalue
    
    taumax_matrix[one_way_index] = all_cor$taumax
    taumax_matrix[back_way_index] = all_cor$taumax
    
    completeness_matrix[one_way_index] = all_cor$completeness
    completeness_matrix[back_way_index] = all_cor$completeness
    
    return(list(cor = cor_matrix, raw = raw_matrix, pvalue = pvalue_matrix, taumax = taumax_matrix, completeness = completeness_matrix, keep = t(!exclude_loc)))
  } else {
    return(list(cor = all_cor))
  }
}

#' Fast kendall tau
#' 
#' Uses the underlying c++ implementation of `ici_kt` to provide a fast version
#' of Kendall-tau correlation.
#' 
#' @param x a numeric vector, matrix, or data frame.
#' @param y NULL (default) or a vector.
#' @param use an optional character string giving a method for computing correlations in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", or "pairwise.complete.obs".
#' @param alternative the type of test
#' @param continuity should a continuity correction be applied
#' @param return_matrix Should the matrices of values be returned, or a long data.frame
#' 
#' @details Although the interface is *mostly* identical to the built-in 
#' [stats::cor()] method, there are some differences. 
#'   
#'  * if only `x` is provided as a matrix or data.frame, the columns must be named.
#'  * if providing both `x` and `y`, it is assumed they are both
#'  single vectors.
#'  * if `NA` values are present, this function does not error, but will either remove them 
#'  or return `NA`, depending on the option.
#'  * "na.or.complete" is not a valid option for `use`.
#'  * A named list with matrices or data.frame is returned, with the `tau` and `pvalue` values.
#' 
#' @return a list of matrices, tau, pvalue, or a data.frame.
#' @export
kt_fast = function(x, y = NULL, use = "everything",
                  alternative = "two.sided",
                  continuity = FALSE, return_matrix = TRUE)
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
                                        which = "kt",
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
      in_split$tau = NA
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
    split_cor = computation$split_fun(split_comparisons, kt_split, x, na_method, do_log_memory)
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
    tau_matrix = cor_matrix
    pvalue_matrix = cor_matrix
    
    one_way_index = cbind(all_cor$s1, all_cor$s2)
    back_way_index = cbind(all_cor$s2, all_cor$s1)
    
    tau_matrix[one_way_index] = all_cor$tau
    tau_matrix[back_way_index] = all_cor$tau
    
    pvalue_matrix[one_way_index] = all_cor$pvalue
    pvalue_matrix[back_way_index] = all_cor$pvalue
    
    return(list(tau = tau_matrix, pvalue = pvalue_matrix, run_time = t_diff))
  } else {
    return(list(tau = all_cor, run_time = t_diff))
  }
  
}


#' pairwise completeness
#' 
#' Calculates the completeness between any two samples using "or", is an
#' entry missing in either X "or" Y.
#' 
#' @param data_matrix samples are columns, features are rows
#' @param global_na globally, what should be treated as NA?
#' @param include_only is there certain comparisons to do?
#' @param return_matrix should the matrix or data.frame be returned?
#' 
#' @seealso [ici_kendalltau()]
#' 
#' @export
#' 
#' @return matrix of degree of completeness
pairwise_completeness = function(data_matrix,
                                global_na = c(NA, Inf, 0),
                                include_only = NULL,
                                return_matrix = TRUE){
  
  data_arg = rlang::caller_arg(data_matrix)
  include_arg = rlang::caller_arg(include_only)
  
  check_if_colnames_null(data_matrix, arg = data_arg)
  data_matrix = transform_to_matrix(data_matrix, arg = data_arg)
  check_if_numeric(data_matrix, arg = data_arg)
  
  
  exclude_loc = setup_missing_matrix(data_matrix, global_na)
  
  n_sample = ncol(exclude_loc)
  
  computation = check_furrr()

  split_comparisons = setup_comparisons(colnames(data_matrix),
                                        include_only = include_only,
                                        diag_good = FALSE,
                                        ncore = computation$ncore,
                                        which = "completeness",
                                        include_arg = include_arg)
  
  
  
  split_missing = computation$split_fun(split_comparisons, complete_split, exclude_loc)
  
  all_missing = purrr::list_rbind(split_missing)
  rownames(all_missing) = NULL
  
  if (return_matrix) {
    completeness_matrix = matrix(0, nrow = ncol(exclude_loc), ncol = ncol(exclude_loc))
    rownames(completeness_matrix) = colnames(completeness_matrix) = colnames(exclude_loc)
    
    completeness_matrix[cbind(all_missing$s1, all_missing$s2)] = all_missing$completeness
    completeness_matrix[cbind(all_missing$s2, all_missing$s1)] = all_missing$completeness
    
    return(completeness_matrix)
  } else {
    return(all_missing)
  }

}

complete_split = function(do_comparisons, exclude_loc) {
  #seq_range = seq(in_range[1], in_range[2])
  #print(seq_range)
  missingness = vector("numeric", nrow(do_comparisons))
  completeness = vector("numeric", nrow(do_comparisons))
  for (irow in seq(1, nrow(do_comparisons))) {
    iloc = do_comparisons[irow, 1]
    jloc = do_comparisons[irow, 2]
    missingness[irow] = missing_either(exclude_loc[, iloc], exclude_loc[, jloc])
  }
  completeness = 1 - (missingness / nrow(exclude_loc))
  do_comparisons$missingness = missingness
  do_comparisons$completeness = completeness
  do_comparisons
}

missing_either = function(in_x, in_y){
  not_in_one = sum(in_x | in_y)
  not_in_one
}



check_icikt_timing = function(exclude_data, tmp_pairwise, perspective, n_todo, ncore){
  t_start = Sys.time()
  for (irow in seq_len(nrow(tmp_pairwise))) {
    iloc = tmp_pairwise[irow, 1]
    jloc = tmp_pairwise[irow, 2]
    tmp_val = ici_kt(exclude_data[, iloc], exclude_data[, jloc], perspective = perspective)
  }
  t_stop = Sys.time()
  t_total = as.numeric(difftime(t_stop, t_start, units = "secs"))
  n_comp = nrow(tmp_pairwise)
  
  t_each = t_total / n_comp
  
  t_theoretical = t_each * n_todo
  t_cores = t_theoretical / ncore
  
  data.frame(which =
               c("n_tested",
                 "n_todo",
                 "time_tested",
                 "time_single",
                 "time_all",
                 "time_across_cores",
                 "time_minutes",
                 "time_hours",
                 "time_days"),
             value = c(n_comp,
                n_todo,
                t_total,
                t_each,
                t_theoretical,
                t_cores,
                t_cores / 60,
                t_cores / (60 * 60),
                t_cores / (60 * 60 * 60)))
  
}

